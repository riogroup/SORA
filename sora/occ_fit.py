from .config import test_attr
from .star import Star
from .ephem import Ephemeris, EphemPlanete, EphemJPL, EphemKernel
from .observer import Observer
from .lightcurve import LightCurve
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


class _PositionDict(dict):
    def __setitem__(self, key, value):
        if key.startswith('_occ_'):
            super().__setitem__(key[5:], value)
        elif key in self.keys():
            status = {'on': True, 'off': False}
            if value in status.keys():
                if type(self[key]) == _PositionDict:
                    for k in self[key].keys():
                        self[key][k] = value
                elif key=='on':
                    super().__setitem__('on', status[value])
            else:
                raise ValueError("Value must be 'on' or 'off' only.")
        else:
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
        if len(self.__observations) > 0:
            for o,l in self.__observations:
                if o == obs and l == lightcurve:
                    raise ValueError('{} observation already defined'.format(obs.name))
        self.__observations.append((obs,lightcurve))
        lightcurve.set_vel(np.absolute(self.vel))
        lightcurve.set_dist(float(self.dist.AU))
        lightcurve.set_diam(float(self.star_diam.AU))

    def fit_ellipse(self):
        # fit ellipse to the points
        return

    def fit_to_shape(self):
        # fit points to a 3D shape model
        return

    @property
    def positions(self):
        position = _PositionDict()
        if hasattr(self,'_position'):
            position = self._position
        if len(self.__observations) == 0:
            raise ValueError('There is no observation defined for this occultation')

        for o,l in self.__observations:
            if not o.name in position.keys():
                position['_occ_'+o.name] = _PositionDict()
            status = False
            if hasattr(l, 'immersion'):
                if not 'immersion' in position[o.name].keys():
                    position[o.name]['_occ_immersion'] = _PositionDict(on=True)
                obs_im = position[o.name]['immersion']
                if 'time' in obs_im.keys() and obs_im['time'] == l.immersion:
                    pass
                else:
                    f1,g1 = positionv(self.star,self.ephem,o,l.immersion)[0:2]
                    fe1,ge1 = positionv(self.star,self.ephem,o,l.immersion-l.immersion_err*u.s)[0:2]
                    fe2,ge2 = positionv(self.star,self.ephem,o,l.immersion+l.immersion_err*u.s)[0:2]
                    obs_im['_occ_time'] =l.immersion
                    #obs_im['_occ_err'] = l.immersion_err
                    obs_im['_occ_value'] = (round(f1,3),round(g1,3))
                    obs_im['_occ_error'] = ((round(fe1,3),round(ge1,3)),(round(fe2,3),round(ge2,3)))
                status = True
            if hasattr(l, 'emersion'):
                if not 'emersion' in position[o.name].keys():
                    position[o.name]['_occ_emersion'] = _PositionDict(on=True)
                obs_em = position[o.name]['emersion']
                if 'time' in obs_em.keys() and obs_em['time'] == l.emersion:
                    pass
                else:
                    f1,g1 = positionv(self.star,self.ephem,o,l.emersion)[0:2]
                    fe1,ge1 = positionv(self.star,self.ephem,o,l.emersion-l.emersion_err*u.s)[0:2]
                    fe2,ge2 = positionv(self.star,self.ephem,o,l.emersion+l.emersion_err*u.s)[0:2]
                    obs_em['_occ_time'] =l.emersion
                    #obs_em['_occ_err'] = l.emersion_err
                    obs_em['_occ_value'] = (round(f1,3),round(g1,3))
                    obs_em['_occ_error'] = ((round(fe1,3),round(ge1,3)),(round(fe2,3),round(ge2,3)))
                status = True
            if not status:
                f1,g1 = positionv(self.star,self.ephem,o,l.initial_time)[0:2]
                position[o.name]['_occ_start_obs'] = (f1,g1)
                f1,g1 = positionv(self.star,self.ephem,o,l.end_time)[0:2]
                position[o.name]['_occ_end_obs'] = (f1,g1)
            if status:
                position[o.name]['_occ_status'] = 'positive'
            else:
                position[o.name]['_occ_status'] = 'negative'
        self._position = position
        return self._position

    def plot_chords(self, positive_color='blue', negative_color='green', error_color='red'):
        # plot chords of the occultation
        positions = self.positions
        for site in positions.keys():
            if positions[site]['status'] == 'negative':
                arr = np.array([positions[site]['start_obs'], positions[site]['end_obs']])
                plt.plot(*arr.T, '--', color=negative_color, linewidth=0.7)
            else:
                n = 0
                if positions[site]['immersion']['on']:
                    arr = np.array([positions[site]['immersion']['error']])
                    plt.plot(*arr.T, color=error_color, linewidth=1.5)
                    n+=1
                if positions[site]['emersion']['on']:
                    arr = np.array([positions[site]['emersion']['error']])
                    plt.plot(*arr.T, color=error_color, linewidth=1.5)
                    n+=1
                if n == 2:
                    arr = np.array([positions[site]['immersion']['value'], positions[site]['emersion']['value']])
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