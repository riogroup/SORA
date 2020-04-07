from .config import test_attr
from .star import Star
from .ephem import Ephemeris, EphemPlanete, EphemJPL, EphemKernel
from .observer import Observer
from .lightcurve import LightCurve
from .extra import ChiSquare
import astropy.units as u
import astropy.constants as const
from astropy.time import Time
from astropy.coordinates import SphericalCosLatDifferential, SkyCoord, SkyOffsetFrame
from astropy.table import Table
from astroquery.vizier import Vizier
import numpy as np
import warnings
import matplotlib.pyplot as plt


def prediction(ephem, time_beg, time_end, mag_lim=None, interv=60, divs=1):
    """ Predicts stellar occultations
        
    Parameters:
    ephem (Ephem): Ephemeris. It must be an Ephemeris object.
    time_beg (Time): Initial time for prediction
    time_beg (Time): Final time for prediction
    mag_lim (int,float): Faintest Gmag for search
    interv (int, float): interval, in seconds, of ephem times for search
    divs (int,float): interal, in deg, for max search of stars
    
    Return:
    occ_params (Table): Astropy Table with the occultation params for each event
    """
    # generate ephemeris
    if type(ephem) != EphemKernel:
        raise TypeError('At the moment prediction only works with EphemKernel')
    print("Generating Ephemeris ...")
    dt = np.arange(0, (time_end-time_beg).sec, interv)*u.s
    t = time_beg + dt
    coord = ephem.get_position(t)

    # define catalogue parameters
    kwds = {}
    kwds['columns'] = ['Source', 'RA_ICRS', 'DE_ICRS']
    kwds['row_limit']=10000000
    kwds['timeout']=600
    if mag_lim:
        kwds['column_filters']={"Gmag":"<{}".format(mag_lim)}
    vquery = Vizier(**kwds)

    # determine suitable divisions for star search
    radius = ephem.radius + const.R_earth
    mindist = np.arcsin(radius/coord[0].distance)
    divisions = []
    n=0
    while True: 
        dif = coord.separation(coord[n]) 
        k = np.where(dif < divs*u.deg)[0] 
        l = np.where(k[1:]-k[:-1] > 1)[0] 
        l = np.append(l,len(k)-1) 
        m = np.where(k[l] - n > 0)[0].min() 
        k = k[l][m] 
        divisions.append([n,k])
        if k == len(coord)-1: 
            break 
        n = k 
        
    print('Ephemeris was split in {} parts for better search of stars'.format(len(divisions)))

    # makes predictions for each division
    occs = []
    for i,vals in enumerate(divisions):
        print('\nSearching occultations in part {}/{}'.format(i+1,len(divisions)))
        nt = t[vals[0]:vals[1]]
        ncoord = coord[vals[0]:vals[1]]
        ra = np.mean([ncoord.ra.min().deg,ncoord.ra.max().deg])
        dec = np.mean([ncoord.dec.min().deg,ncoord.dec.max().deg])
        width = ncoord.ra.max() - ncoord.ra.min() + 2*mindist
        height = ncoord.dec.max() - ncoord.dec.min() + 2*mindist
        pos_search = SkyCoord(ra*u.deg, dec*u.deg)
        
        print('Downloading stars ...')
        catalogue = vquery.query_region(pos_search, width=width, height=height, catalog='I/345/gaia2')
        print('Identifying occultations ...')
        if len(catalogue) == 0:
            continue
        catalogue = catalogue[0]
        stars = SkyCoord(catalogue['RA_ICRS'], catalogue['DE_ICRS'])
        idx, d2d, d3d = stars.match_to_catalog_sky(ncoord)
        
        dist = np.arcsin(radius/ncoord[idx].distance)
        k = np.where(d2d < dist)[0]
        for ev in k:
            star = Star(code=catalogue['Source'][ev], log=False)
            pars = [star.code, star.geocentric(nt[idx][ev]), star.mag['G']]
            pars = np.hstack((pars, occ_params(star, ephem, nt[idx][ev])))
            occs.append(pars)

    if not occs:
        warnings.warn('No stellar occultation was found')
        return Table(names=['time', 'coord', 'ca', 'pa', 'vel', 'G', 'G*', 'dist'])
    # create astropy table with the params
    occs2 = np.transpose(occs)
    time = Time(occs2[3])
    k = np.argsort(time)
    source = occs2[0][k]
    coord = [i.to_string('hmsdms',precision=5, sep=' ') for i in occs2[1][k]]
    mags = ['{:6.3f}'.format(i) for i in occs2[2][k]]
    magstt = ['{:6.3f}'.format(occs2[2][i] + 2.5*np.log10(np.absolute(occs2[6][i].value)/20.0)) for i in k]
    time = [i.iso for i in time[k]]
    ca = ['{:5.3f}'.format(i.value) for i in occs2[4][k]]
    pa = ['{:6.2f}'.format(i.value) for i in occs2[5][k]]
    vel = ['{:-6.2f}'.format(i.value) for i in occs2[6][k]]
    dist = ['{:7.3f}'.format(i.value) for i in occs2[7][k]]
    t = Table([time, coord, ca, pa, vel, mags, magstt, dist],
               names=['time', 'coord', 'ca', 'pa', 'vel', 'G', 'G*', 'dist'])
    return t


def positionv(star,ephem,observer,time):
    """ Calculates the position and velocity of the occultation shadow relative to the observer.
        
    Parameters:
    star (Star): The coordinate of the star in the same frame as the ephemeris.
    It must be a Star object.
    ephem (Ephem): Ephemeris. It must be an Ephemeris object.
    observer (Observer): The Observer object to be added.
    time (Time): Instant to calculate position and velocity
    
    Return:
    f, g (float): The orthographic projection of the shadow relative to the observer
    """
    if type(star) != Star:
        raise ValueError('star must be a Star object')
    if type(ephem) not in [Ephemeris, EphemPlanete, EphemJPL, EphemKernel]:
        raise ValueError('ephem must be an Ephemeris object')
    if type(observer) != Observer:
        raise ValueError('observer must be an Observer object')
        
    coord = star.geocentric(time)
    dt = 0.1*u.s
    
    if type(ephem) == EphemPlanete:
        ephem.fit_d2_ksi_eta(coord, log=False)
    ksio1, etao1 = observer.get_ksi_eta(time=time, star=coord)
    ksie1, etae1 = ephem.get_ksi_eta(time=time, star=coord)
    
    f = ksio1+ksie1
    g = etao1+etae1
    
    ksio2, etao2 = observer.get_ksi_eta(time=time+dt, star=coord)
    ksie2, etae2 = ephem.get_ksi_eta(time=time+dt, star=coord)
    
    nf = ksio2+ksie2
    ng = etao2+etae2
    
    vf = (nf-f)/0.1
    vg = (ng-g)/0.1

    return f, g, vf, vg


def occ_params(star, ephem, time):
    """ Calculates the parameters of the occultation, as instant, CA, PA.
        
    Parameters:
    star (Star): The coordinate of the star in the same frame as the ephemeris.
    It must be a Star object.
    ephem (Ephem): Ephemeris. It must be an Ephemeris object.
    
    Return:
    instant of CA (Time): Instant of Closest Approach
    CA (arcsec): Distance of Closest Approach
    PA (deg): Position Angle at Closest Approach
    """
    
    delta_t = 0.05
    
    if type(star) != Star:
        raise ValueError('star must be a Star object')
    if type(ephem) not in [Ephemeris, EphemKernel, EphemJPL, EphemPlanete]:
        raise ValueError('ephem must be a Ephemeris object')
        
    tt = time + np.arange(-600, 600, delta_t)*u.s
    coord = star.geocentric(tt[0])
    if type(ephem) == EphemPlanete:
        ephem.fit_d2_ksi_eta(coord, log=False)
    ksi, eta = ephem.get_ksi_eta(tt, coord)
    dd = np.sqrt(ksi*ksi+eta*eta)
    min = np.argmin(dd)
    
    if type(ephem) == EphemPlanete:
        dist = ephem.ephem[int(len(ephem.time)/2)].distance
    else:
        dist = ephem.get_position(time).distance
    
    ca = np.arcsin(dd[min]*u.km/dist).to(u.arcsec)
    
    pa = (np.arctan2(-ksi[min],-eta[min])*u.rad).to(u.deg)
    if pa < 0*u.deg:
        pa = pa + 360*u.deg
    
    dksi = ksi[min+1]-ksi[min]
    deta = eta[min+1]-eta[min]
    vel = np.sqrt(dksi**2 + deta**2)/delta_t
    vel = -vel*np.sign(dksi)*(u.km/u.s)
    
    return tt[min], ca, pa, vel, dist.to(u.AU)


def fit_ellipse(*args, **kwargs):
    params_needed = ['center_f', 'center_g', 'equatorial_radius', 'oblateness', 'pos_angle']
    if not all([param in kwargs for param in params_needed]):
        raise ValueError('Input conditions not satisfied. Please refer to the tutorial.')
    center_f = kwargs['center_f']
    dcenter_f = 0.0
    if 'dcenter_f' in kwargs:
        dcenter_f = kwargs['dcenter_f']
    center_g = kwargs['center_g']
    dcenter_g = 0.0
    if 'dcenter_g' in kwargs:
        dcenter_g = kwargs['dcenter_g']
    equatorial_radius = kwargs['equatorial_radius']
    dequatorial_radius = 0.0
    if 'dequatorial_radius' in kwargs:
        dequatorial_radius = kwargs['dequatorial_radius']
    oblateness = kwargs['oblateness']
    doblateness = 0.0
    if 'doblateness' in kwargs:
        doblateness = kwargs['doblateness']
    pos_angle = kwargs['pos_angle']
    dpos_angle = 0.0
    if 'dpos_angle' in kwargs:
        dpos_angle = kwargs['dpos_angle']
    loop = 10000000
    if 'loop' in kwargs:
        loop = kwargs['loop']
    number_chi = 10000
    if 'number_chi' in kwargs:
        number_chi = kwargs['number_chi']
    log = False
    if 'log' in kwargs:
        log = kwargs['log']

    values = []
    for occ in args:
        if type(occ) != Occultation:
            raise TypeError('Given argument must be an Occultation object.')
        pos = occ.positions
        for site in pos.keys():
            pos_obs = pos[site]
            for lc in pos_obs.keys():
                pos_lc = pos_obs[lc]
                if type(pos_lc) != _PositionDict:
                    continue
                if pos_lc['status'] == 'positive':
                    if pos_lc['immersion']['on']:
                        f,g = pos_lc['immersion']['value']
                        err = np.array(pos_lc['immersion']['error'])
                        erro = np.linalg.norm(err[0]-err[1])/2.0
                        values.append([f,g,erro])
                    if pos_lc['emersion']['on']:
                        f,g = pos_lc['emersion']['value']
                        err = np.array(pos_lc['emersion']['error'])
                        erro = np.linalg.norm(err[0]-err[1])/2.0
                        values.append([f,g,erro])

    controle_f0 = Time.now()
    f0_chi = np.array([])
    g0_chi = np.array([])
    a_chi  = np.array([])
    obla_chi  = np.array([])
    posang_chi  = np.array([])
    chi2_best = np.array([])

    while (len(f0_chi) < number_chi):
        chi2 = np.zeros(loop)
        f0 = center_f + dcenter_f*(2*np.random.random(loop) - 1)
        g0 = center_g + dcenter_g*(2*np.random.random(loop) - 1)
        a  = equatorial_radius + dequatorial_radius*(2*np.random.random(loop) - 1)
        obla  = oblateness + doblateness*(2*np.random.random(loop) - 1)
        obla[obla<0],obla[obla>1] = 0, 1
        phi_deg  = pos_angle + dpos_angle*(2*np.random.random(loop) - 1)
        controle_f1 = Time.now()

        for fi, gi, si in values:
            b = a - a*obla
            phi = phi_deg*(np.pi/180.0)
            dfi = fi-f0
            dgi = gi-g0
            r = np.sqrt(dfi**2 + dgi**2)
            theta = np.arctan2(dgi,dfi)
            ang = theta+phi
            r_model = (a*b)/np.sqrt((a*np.sin(ang))**2 + (b*np.cos(ang))**2)
            f_model = f0 + r_model*np.cos(theta)
            g_model = g0 + r_model*np.sin(theta)
            chi2 += ((fi - f_model)**2 + (gi - g_model)**2)/(si**2)

        controle_f2 = Time.now()
        if 'dchi_min' in kwargs:
            region = np.where(chi2 < chi2.min() + kwargs['dchi_min'])[0]
        else:
            region = np.arange(len(chi2))
        chi2_best = np.append(chi2_best,chi2[region])
        if log:
            print('Elapsed time: {:.3f} seconds.'.format((controle_f2 - controle_f1).sec))
            print(len(chi2[region]),len(chi2_best)) #TEST
        f0_chi = np.append(f0_chi,f0[region])
        g0_chi = np.append(g0_chi,g0[region])
        a_chi  = np.append(a_chi,a[region])
        obla_chi  = np.append(obla_chi,obla[region])
        posang_chi  = np.append(posang_chi,phi_deg[region])

    chisquare = ChiSquare(chi2_best, len(values), center_f=f0_chi, center_g=g0_chi, equatorial_radius=a_chi,
                          oblateness=obla_chi, position_angle=posang_chi)
    controle_f4 = Time.now()
    if log:
        print('Total elapsed time: {:.3f} seconds.'.format((controle_f4 - controle_f0).sec))
    return chisquare


class _PositionDict(dict):
    def __setitem__(self, key, value):
        status = {'on': True, 'off': False}
        n = 0
        if key.startswith('_occ_'):
            super().__setitem__(key[5:], value)
            n = 1
        elif key in self.keys():
            n = 1
            if value not in status.keys():
                raise ValueError("Value must be 'on' or 'off' only.")
            if type(self[key]) == _PositionDict:
                for k in self[key].keys():
                    self[key][k] = value
            elif key=='on':
                super().__setitem__('on', status[value])
        else:
            if value not in status.keys():
                raise ValueError("Value must be 'on' or 'off' only.")
            for key1 in self.keys():
                if type(self[key1]) == _PositionDict:
                    if key in self[key1].keys():
                        n = 1
                        self[key1][key] = value
        if n==0:
            raise KeyError('Key "{}" does not exist'.format(key))

    def __str__(self):
        out = '\n' + '\n'.join(['{}: {}'.format(key, self[key]) for key in self.keys()])
        return out.replace('\n', '\n  ')


### Object for occultation
class Occultation():
    '''
    Docstring
    Do the reduction of the occultation
    '''
    def __init__(self, star, ephem, time):
        """ Instantiate Occultation object.
        
        Parameters:
        star (Star):The coordinate of the star in the same frame as the ephemeris.
        It must be a Star object.
        ephem (Ephem):Ephemeris. It must be an Ephemeris object.

        """
        if type(star) != Star:
            raise ValueError('star must be a Star object')
        if type(ephem) not in [Ephemeris, EphemKernel, EphemJPL]:
            raise ValueError('ephem must be a Ephemeris object')
        self.star = star
        self.ephem = ephem
        
        tt, ca, pa, vel, dist = occ_params(star,ephem, time)
        self.ca = ca   # Closest Approach distance
        self.pa = pa   # Position Angle at CA
        self.vel = vel  # Shadow velocity at CA
        self.dist = dist  # object distance at CA
        self.tca = tt   # Instant of CA
        self.star_diam = self.star.apparent_diameter(self.dist, log=False)
        
        self.__observations = []
        self._position = _PositionDict()
    
    def add_observation(self, obs, lightcurve):
        """ Add Observers to the Occultation object.
        
        Parameters:
        obs (Observer):The Observer object to be added.
        status (string): it can be "positive", "negative", "visual" or "undefined"

        """
        if type(obs) != Observer:
            raise ValueError('obs must be an Observer object')
        if type(lightcurve) != LightCurve:
            raise ValueError('lightcurve must be a LightCurve object')
        for o,l in self.__observations:
            if l.name == lightcurve.name:
                raise ValueError('{} LightCurve already associated to {} Observer'.format(lightcurve.name, o.name))
        self.__observations.append((obs,lightcurve))
        lightcurve.set_vel(np.absolute(self.vel))
        lightcurve.set_dist(float(self.dist.AU))
        lightcurve.set_diam(float(self.star_diam.AU))

    def remove_observation(self, key, key_lc=None):
        rm_list = np.array([])
        obs = []
        lcs = []
        same_key = False
        if key_lc is None:
            same_key = True
            key_lc = key
        ko = []
        kl = []
        for i, val in enumerate(self.__observations):
            if val[0].name == key:
                ko.append(i)
            if val[1].name == key_lc:
                kl.append(i)
        if not same_key:
            k  = np.where(np.array(ko) == np.array(kl))[0]
            rm_list = np.hstack((rm_list, k))
        else:
            rm_list = np.hstack((rm_list, np.array(kl)))
            if len(kl) > 0 and len(ko) > 0 and kl[0] not in ko:
                raise ValueError("Observation could not univocally be identified, please give parameters for Observer and LightCurve")
            rm_list = np.hstack((rm_list, np.array(ko)))
        rm_list = np.unique(rm_list)
        if len(rm_list) == 0:
            raise ValueError('No observer "{}" and/or lightcurve "{}" was found'.format(key,key_lc))
        list = np.arange(len(self.__observations)).tolist()
        for i in rm_list:
            list.remove(i)
        self.__observations = [self.__observations[item] for item in list]

    def observations(self):
        """ Print all the observations added to the Occultation object
        Pair (Observer, LightCurve)
        """
        for o,l in self.__observations:
            print('Observer= {}, LC: {}'.format(o.name, l.name))

    def fit_ellipse(self, **kwargs):
        # fit ellipse to the points
        chisquare = fit_ellipse(self, **kwargs)
        return chisquare

    def fit_to_shape(self):
        # fit points to a 3D shape model
        return

    @property
    def positions(self):
        position = self._position
        if len(self.__observations) == 0:
            raise ValueError('There is no observation defined for this occultation')

        pair = []
        for o,l in self.__observations:
            pair.append((o.name, l.name))

            coord = [o.lon,o.lat,o.height]
            if o.name not in position.keys():
                position['_occ_'+o.name] = _PositionDict(lon=o.lon, lat=o.lat, height=o.height)
                position[o.name]['_occ_lon'] = o.lon
                position[o.name]['_occ_lat'] = o.lat
                position[o.name]['_occ_height'] = o.height
            pos_obs = position[o.name]
            coord2 = [pos_obs['lon'],pos_obs['lat'],pos_obs['height']]
            if o.lon != pos_obs['lon']:
                position[o.name]['_occ_lon'] = o.lon
            if o.lat != pos_obs['lat']:
                position[o.name]['_occ_lat'] = o.lat
            if o.height != pos_obs['height']:
                position[o.name]['_occ_height'] = o.height
            samecoord = (coord == coord2)

            if l.name not in pos_obs.keys():
                pos_obs['_occ_'+l.name] = _PositionDict()
            pos_lc = pos_obs[l.name]

            pos_lc['_occ_status'] = 'negative'
            if hasattr(l, 'immersion') or hasattr(l, 'emersion'):
                pos_lc['_occ_status'] = 'positive'

            if hasattr(l, 'immersion'):
                if 'immersion' not in pos_lc.keys():
                    pos_lc['_occ_immersion'] = _PositionDict(on=True)
                obs_im = pos_lc['immersion']
                do_err = False
                if samecoord and 'time' in obs_im.keys() and obs_im['time'] == l.immersion:
                    pass
                else:
                    do_err = True
                    f1,g1 = positionv(self.star,self.ephem,o,l.immersion)[0:2]
                    obs_im['_occ_time'] =l.immersion
                    obs_im['_occ_value'] = (round(f1,3),round(g1,3))
                if not do_err and 'time_err' in obs_im.keys() and obs_im['time_err'] == l.immersion_err:
                    pass
                else:
                    fe1,ge1 = positionv(self.star,self.ephem,o,l.immersion-l.immersion_err*u.s)[0:2]
                    fe2,ge2 = positionv(self.star,self.ephem,o,l.immersion+l.immersion_err*u.s)[0:2]
                    obs_im['_occ_time_err'] = l.immersion_err
                    obs_im['_occ_error'] = ((round(fe1,3),round(ge1,3)),(round(fe2,3),round(ge2,3)))

            if hasattr(l, 'emersion'):
                if 'emersion' not in pos_lc.keys():
                    pos_lc['_occ_emersion'] = _PositionDict(on=True)
                obs_em = pos_lc['emersion']
                do_err = False
                if samecoord and 'time' in obs_em.keys() and obs_em['time'] == l.emersion:
                    pass
                else:
                    do_err = True
                    f1,g1 = positionv(self.star,self.ephem,o,l.emersion)[0:2]
                    obs_em['_occ_time'] =l.emersion
                    obs_em['_occ_value'] = (round(f1,3),round(g1,3))
                if not do_err and 'time_err' in obs_em.keys() and obs_em['time_err'] == l.emersion_err:
                    pass
                else:
                    fe1,ge1 = positionv(self.star,self.ephem,o,l.emersion-l.emersion_err*u.s)[0:2]
                    fe2,ge2 = positionv(self.star,self.ephem,o,l.emersion+l.emersion_err*u.s)[0:2]
                    obs_em['_occ_time_err'] = l.emersion_err
                    obs_em['_occ_error'] = ((round(fe1,3),round(ge1,3)),(round(fe2,3),round(ge2,3)))

            if pos_lc['status'] == 'negative':
                if 'start_obs' not in pos_lc.keys():
                    pos_lc['_occ_start_obs'] = _PositionDict(on=True)
                obs_start = pos_lc['start_obs']
                if samecoord and 'time' in obs_start.keys() and obs_start['time'] == l.initial_time:
                    pass
                else:
                    f1,g1 = positionv(self.star,self.ephem,o,l.initial_time)[0:2]
                    obs_start['_occ_time'] = l.initial_time
                    obs_start['_occ_value'] = (f1,g1)
                if 'end_obs' not in pos_lc.keys():
                    pos_lc['_occ_end_obs'] = _PositionDict(on=True)
                obs_end = pos_lc['end_obs']
                if samecoord and 'time' in obs_end.keys() and obs_end['time'] == l.end_time:
                    pass
                else:
                    f1,g1 = positionv(self.star,self.ephem,o,l.end_time)[0:2]
                    obs_end['_occ_time'] = l.end_time
                    obs_end['_occ_value'] = (f1,g1)

        for key in list(position):
            for key_lc in list(position[key]):
                if type(key_lc) == _PositionDict and (key,key_lc) not in pair:
                    del position[key][key_lc]
            if len(position[key]) == 0:
                del position[key]

        return self._position

    @positions.setter
    def positions(self, value):
        if not hasattr(self, '_position'):
            pos = self.positions
        if value not in ['on','off']:
            raise ValueError("Value must be 'on' or 'off' only.")
        for key in self._position.keys():
            self._position[key] = value

    def plot_chords(self, all_chords=True, positive_color='blue', negative_color='green', error_color='red'):
        # plot chords of the occultation
        positions = self.positions
        for site in positions.keys():
            pos_obs = positions[site]
            for lc in pos_obs.keys():
                pos_lc = pos_obs[lc]
                if type(pos_lc) != _PositionDict:
                    continue
                if pos_lc['status'] == 'negative':
                    arr = np.array([pos_lc['start_obs']['value'], pos_lc['end_obs']['value']])
                    plt.plot(*arr.T, '--', color=negative_color, linewidth=0.7)
                else:
                    n = 0
                    if pos_lc['immersion']['on'] or all_chords:
                        arr = np.array([pos_lc['immersion']['error']])
                        plt.plot(*arr.T, color=error_color, linewidth=1.5)
                        n+=1
                    if pos_lc['emersion']['on'] or all_chords:
                        arr = np.array([pos_lc['emersion']['error']])
                        plt.plot(*arr.T, color=error_color, linewidth=1.5)
                        n+=1
                    if n == 2:
                        arr = np.array([pos_lc['immersion']['value'], pos_lc['emersion']['value']])
                        plt.plot(*arr.T, color=positive_color, linewidth=0.7)
        plt.axis('equal')
    
    def plot_occ_map(self):
        # plot occultation map
        return

    def __str__(self):
        """String representation of the Star class
        """
        out = 'Stellar occultation of star Gaia-DR2 {} by {}.\n\n'.format(self.star.code, self.ephem.name)
        out += 'Geocentric Closest Approach: {:.3f}\n'.format(self.ca)
        out += 'Instant of CA: {}\n'.format(self.tca.iso)
        out += 'Position Angle: {:.2f}\n'.format(self.pa)
        out += 'Geocentric shadow velocity: {:.2f}\n\n'.format(self.vel)

        out += self.star.__str__() + '\n'
        out += self.ephem.__str__() + '\n'

        return out

        self.__count = 0
        self.__out1 = ''
        n = 0
        out1 = ''
        
        n += len(self.obs_positive)
        if len(self.obs_positive) > 0:
            out1 += '{} positive observations\n'.format(len(self.obs_positive))
            for i in self.obs_positive:
                out1 += i.__str__() + '\n'
                out1 += '\n'
            out1 += '\n'
        
        n += len(self.obs_negative)
        if len(self.obs_negative) > 0:
            out1 += '{} negative observations\n'.format(len(self.obs_negative))
            for i in self.obs_negative:
                out1 += i.__str__() + '\n'
                out1 += '\n'
            out1 += '\n'
        
        n += len(self.obs_visual)
        if len(self.obs_visual) > 0:
            out1 += '{} visual observations\n'.format(len(self.obs_visual))
            for i in self.obs_visual:
                out1 += i.__str__() + '\n'
                out1 += '\n'
            out1 += '\n'
        
        n += len(self.obs_undefined)
        if len(self.obs_undefined) > 0:
            out1 += '{} without status observations\n'.format(len(self.obs_undefined))
            for i in self.obs_undefined:
                out1 += i.__str__() + '\n'
                out1 += '\n'
            out1 += '\n'
        out1 += '\b\b'
        
        if n == 0:
            out += 'No observations reported'
        else:
            out += '{} observations reported\n\n'.format(n)
            out += out1
        
        return out