from .config import test_attr
from .star import Star
from .ephem import EphemPlanete, EphemJPL, EphemKernel
from .observer import Observer
from .lightcurve import LightCurve
from .prediction import occ_params, Prediction
from .extra import ChiSquare
import astropy.units as u
from astropy.coordinates import SkyCoord, SkyOffsetFrame
from astropy.time import Time
import numpy as np
import warnings
import matplotlib.pyplot as plt


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
    if type(ephem) not in [EphemPlanete, EphemJPL, EphemKernel]:
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


def fit_ellipse(*args, **kwargs):
    """ Fits an ellipse to given occultation

    Parameters:
        required params:
        Each occultation is added as the first arguments directly.
        center_f (int,float): The coordinate in f of the ellipse.
        center_g (int,float): The coordinate in g of the ellipse.
        equatorial_radius (int,float): The Equatorial radius of the ellipse.
        oblateness (int,float): The oblateness of the ellipse.
        pos_angle (int,float): The pole position angle of the ellipse.

        Params for the interval of search, if not given, default is set to zero.
        Search between (value - dvalue) and (value + dvalue):
        dcenter_f (int,float)
        dcenter_g (int,float)
        dequatorial_radius (int,float)
        doblateness (int,float)
        dpos_angle (int,float)

        loop: The number of ellipsis to attempt fitting. Default: 10,000,000
        dchi_min: If given, it will only save ellipsis which chi square are
            smaller than chi_min + dchi_min.
        number_chi: if dchi_min is given, the procedure is repeated until
            number_chi is reached. Default: 10,000
        log: If True, it prints information while fitting. Default: False.

    Return:
        chisquare: A ChiSquare object with all parameters.

    Example:
        fit_ellipse(occ1, **kwargs) to fit the ellipse to the chords of occ1 Occultation object
        fit_ellipse(occ1, occ2, **kwargs) to fit the ellipse to the chords
            of occ1 and occ2 Occultation objects together
    """
    params_needed = ['center_f', 'center_g', 'equatorial_radius', 'oblateness', 'pos_angle']
    if not all([param in kwargs for param in params_needed]):
        raise ValueError('Input conditions not satisfied. Please refer to the tutorial.')
    center_f = kwargs['center_f']
    dcenter_f = kwargs.get('dcenter_f', 0.0)
    center_g = kwargs['center_g']
    dcenter_g = kwargs.get('dcenter_g', 0.0)
    equatorial_radius = kwargs['equatorial_radius']
    dequatorial_radius = kwargs.get('dequatorial_radius', 0.0)
    oblateness = kwargs['oblateness']
    doblateness = kwargs.get('doblateness', 0.0)
    pos_angle = kwargs['pos_angle']
    dpos_angle = kwargs.get('dpos_angle', 0.0)
    loop = kwargs.get('loop', 10000000)
    number_chi = kwargs.get('number_chi', 10000)
    log = kwargs.get('log', False)

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

    onesigma = chisquare.get_nsigma(sigma=1)
    for occ in args:
        if type(occ) == Occultation:
            occ.fitted_params = {i:onesigma[i] for i in ['equatorial_radius', 'center_f', 'center_g', 'oblateness', 'position_angle']}
    a = occ.fitted_params['equatorial_radius'][0]
    f0 = occ.fitted_params['center_f'][0]
    g0 = occ.fitted_params['center_g'][0]
    obla = occ.fitted_params['oblateness'][0]
    phi_deg = occ.fitted_params['position_angle'][0]
    radial_dispersion = np.array([])
    error_bar = np.array([])
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
        radial_dispersion = np.append(radial_dispersion,r - r_model)
        error_bar = np.append(error_bar,si)
    occ.chi2_params = {'radial_dispersion': radial_dispersion.std(ddof=1)}
    occ.chi2_params['mean_error'] = [error_bar.mean(),error_bar.std()]
    occ.chi2_params['chi2_min'] = chisquare.get_nsigma()['chi2_min']
    occ.chi2_params['nparam'] = chisquare.nparam
    occ.chi2_params['npts'] = chisquare.npts
    return chisquare


class _PositionDict(dict):
    """ This is a modified Dictionary object to allow for switch on/off of data points,
    it also avoids for user to change data.
    """
    def __setitem__(self, key, value):
        """Redefine how to set a value to a key in the dictionary.
        it only sets a value if the key starts with '_occ_'.
        Otherwise, it only allows for the user to provide 'on' or 'off'
        which is passed only to change the 'on' keyword.
        """
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
            time (str, Time): The approximate time of the occultation (~10min)
                to calculate occultation params.
        """
        if type(star) != Star:
            raise ValueError('star must be a Star object')
        if type(ephem) not in [EphemPlanete, EphemKernel, EphemJPL]:
            raise ValueError('ephem must be a Ephemeris object')
        self.star = star
        self.ephem = ephem
        
        tca, ca, pa, vel, dist = occ_params(star,ephem, time)
        self.ca = ca   # Closest Approach distance
        self.pa = pa   # Position Angle at CA
        self.vel = vel  # Shadow velocity at CA
        self.dist = dist  # object distance at CA
        self.tca = tca   # Instant of CA
        self.star_diam = self.star.apparent_diameter(self.dist, log=False)

        meta = {'name': self.ephem.name, 'radius': self.ephem.radius.to(u.km).value,
            'error_ra': self.ephem.error_ra.to(u.mas).value, 'error_dec': self.ephem.error_dec.to(u.mas).value}
        self.predict = Prediction(time=[tca], coord_star=[self.star.geocentric(tca)],
            coord_obj=[self.ephem.get_position(tca)], ca=[ca.value], pa=[pa.value], vel=[vel.value],
            dist=[dist.value], mag=[self.star.mag['G']], source=[self.star.code], meta=meta)
        
        self.__observations = []
        self._position = _PositionDict()
    
    def add_observation(self, obs, lightcurve):
        """ Add observations to the Occultation object.

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
        lightcurve.set_diam(float(self.star_diam.km))
        try:
            lightcurve.calc_magnitude_drop(mag_star=self.star.mag['G'],mag_obj=self.ephem.apparent_magnitude(self.tca))
        except:
            lightcurve.bottom_flux = 0.0
            warnings.warn('Magnitude drop was not calculated. Using bottom flux as 0.0 instead.')


    def remove_observation(self, key, key_lc=None):
        """ Removes observation from the Occultation object.

        Parameters:
            key: The name given to Observer or LightCurve to remove from the list.
            keylc: In the case where repeated names are present for different observations,
                keylc must be given for the name of the LightCurve
                and key will be used for the name of the Observer.
        """
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
        """ Fits an ellipse to the chords of this occultation

        Parameters:
            required params:
            Each occultation is added as the first arguments directly.
            center_f (int,float): The coordinate in f of the ellipse.
            center_g (int,float): The coordinate in g of the ellipse.
            equatorial_radius (int,float): The Equatorial radius of the ellipse.
            oblateness (int,float): The oblateness of the ellipse.
            pos_angle (int,float): The pole position angle of the ellipse.

            Params for the interval of search, if not given, default is set to zero.
            Search between (value - dvalue) and (value + dvalue):
            dcenter_f (int,float)
            dcenter_g (int,float)
            dequatorial_radius (int,float)
            doblateness (int,float)
            dpos_angle (int,float)

            loop: The number of ellipsis to attempt fitting. Default: 10,000,000
            dchi_min: If given, it will only save ellipsis which chi square are
                smaller than chi_min + dchi_min.
            number_chi: if dchi_min is given, the procedure is repeated until
                number_chi is reached. Default: 10,000
            log: If True, it prints information while fitting. Default: False.

        Return:
            chisquare: A ChiSquare object with all parameters.
        """
        chisquare = fit_ellipse(self, **kwargs)
        return chisquare

    @property
    def positions(self):
        """ Calculates the positions and velocities for all chords.
        Saves it into an _PositionDict object.
        """
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
                    f1,g1,vf1,vg1 = positionv(self.star,self.ephem,o,l.immersion)
                    obs_im['_occ_time'] =l.immersion
                    obs_im['_occ_value'] = (round(f1,3),round(g1,3))
                    obs_im['_occ_vel'] = (round(vf1,3),round(vg1,3))
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
                    f1,g1,vf1,vg1 = positionv(self.star,self.ephem,o,l.emersion)
                    obs_em['_occ_time'] =l.emersion
                    obs_em['_occ_value'] = (round(f1,3),round(g1,3))
                    obs_em['_occ_vel'] = (round(vf1,3),round(vg1,3))
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
        """ If the users tries to set a value to position, it must be 'on' or 'off',
        and it will be assigned to all chords.
        """
        if not hasattr(self, '_position'):
            pos = self.positions
        if value not in ['on','off']:
            raise ValueError("Value must be 'on' or 'off' only.")
        for key in self._position.keys():
            self._position[key] = value

    def check_velocities(self):
        """ Print the current velocity used by the LightCurves and the Radial velocity.
        """
        if hasattr(self, 'fitted_params'):
            center = np.array([self.fitted_params['center_f'][0], self.fitted_params['center_g'][0]])
        else:
            center = np.array([0,0])
        for o,l in self.__observations:
            vals = self.positions[o.name][l.name]
            if all([i not in vals.keys() for i in ['immersion', 'emersion']]):
                continue
            print('{} - Velocity used: {:.3f}'.format(l.name, l.vel))
            if 'immersion' in vals.keys():
                im = vals['immersion']
                delta = np.array(im['value']) - center
                print('    Immersion Radial Velocity: {:.3f}'.
                      format(np.abs(np.dot(np.array(im['vel']), delta)/np.linalg.norm(delta))))
            if 'emersion' in vals.keys():
                em = vals['emersion']
                delta = np.array(em['value']) - center
                print('    Emersion Radial Velocity: {:.3f}'.
                      format(np.abs(np.dot(np.array(em['vel']), delta)/np.linalg.norm(delta))))

    def new_astrometric_position(self, time=None, offset=None, error=None, log=True):
        """ Calculates the new astrometric position for the object given fitted params

        Parameters:
            time: Time to which calculate the position. If not given, it uses the instant at C/A.
            offset (list): Offset to apply to the position. If not given, it uses the params fitted from ellipse.
                Offsets must be a list of 3 values being [X, Y, 'unit'], where 'unit' must be an angular or distance unit.
                Angular units must be in dacosdec, ddec: Ex: [30.6, 20, 'mas'], or [-15, 2, 'arcsec']
                Distance units must be in X and Y: Ex: [100, -200, 'km'], [0.001, 0.002, 'AU']
            error (list): Error bar of the given offset. If not given, it uses the 1-sigma fitted from ellipse.
                Error must be a list of 3 values being [dX, dY, 'unit'], similar to offset.
                It does not need to be in the same unit as offset.
            log (bool): If true, it prints text, else it returns text.
        """
        if time is not None:
            time = Time(time)
        else:
            time = self.tca

        if offset is not None:
            off_ra = offset[0]*u.Unit(offset[2])
            off_dec = offset[1]*u.Unit(offset[2])
            try:
                teste = off_ra.to(u.km)
                dist = True
            except:
                try:
                    teste = off_ra.to(u.arcsec)
                    dist = False
                except:
                    raise ValueError('Offset unit must be a distance or angular value.')
        elif hasattr(self, 'fitted_params'):
            off_ra = self.fitted_params['center_f'][0]*u.km
            off_dec = self.fitted_params['center_g'][0]*u.km
            dist = True
        else:
            warnings.warn('No offset given or found. Using 0.0 instead.')
            off_ra = 0.0*u.mas
            off_dec = 0.0*u.mas
            dist = False

        if error is not None:
            e_off_ra = error[0]*u.Unit(error[2])
            e_off_dec = error[1]*u.Unit(error[2])
            try:
                teste = e_off_ra.to(u.km)
                e_dist = True
            except:
                try:
                    teste = e_off_ra.to(u.arcsec)
                    e_dist = False
                except:
                    raise ValueError('Error unit must be a distance or angular value.')
        elif hasattr(self, 'fitted_params'):
            e_off_ra = self.fitted_params['center_f'][1]*u.km
            e_off_dec = self.fitted_params['center_g'][1]*u.km
            e_dist = True
        else:
            warnings.warn('No error given or found. Using 0.0 instead.')
            e_off_ra = 0.0*u.mas
            e_off_dec = 0.0*u.mas
            e_dist = False

        coord_geo = self.ephem.get_position(time)
        distance = coord_geo.distance.to(u.km)
        coord_frame = SkyOffsetFrame(origin=coord_geo)
        if dist:
            off_ra = np.arctan2(off_ra, distance)
            off_dec = np.arctan2(off_dec, distance)
        if e_dist:
            e_off_ra = np.arctan2(e_off_ra, distance)
            e_off_dec = np.arctan2(e_off_dec, distance)
        new_pos = SkyCoord(lon=off_ra, lat=off_dec, frame=coord_frame)
        new_pos = new_pos.icrs

        error_star = self.star.error_at(self.tca)
        error_ra = error_star[0] + e_off_ra
        error_dec = error_star[1] + e_off_dec

        out = 'Ephemeris offset (km): X = {:.1f} +/- {:.1f}; Y = {:.1f} +/- {:.1f}\n'.format(
              distance*np.sin(off_ra.to(u.mas)).value, distance*np.sin(e_off_ra.to(u.mas)).value,
              distance*np.sin(off_dec.to(u.mas)).value, distance*np.sin(e_off_dec.to(u.mas)).value)
        out += 'Ephemeris offset (mas): da_cos_dec = {:.3f} +/- {:.3f}; d_dec = {:.3f} +/- {:.3f}\n'.format(
              off_ra.to(u.mas).value, e_off_ra.to(u.mas).value, off_dec.to(u.mas).value, e_off_dec.to(u.mas).value)
        out += '\nAstrometric object position at time {}\n'.format(time.iso)
        out += 'RA = {} +/- {:.3f} mas; DEC = {} +/- {:.3f} mas'.format(
               new_pos.ra.to_string(u.hourangle, precision=5, sep=' '), error_ra.to(u.mas).value,
               new_pos.dec.to_string(u.deg, precision=4, sep=' '), error_dec.to(u.mas).value)

        if log:
            print(out)
        else:
            return out

    def plot_chords(self, all_chords=True, positive_color='blue', negative_color='green', error_color='red'):
        """Plots the chords of the occultation

        Parameters:
            all_chords (bool): if True, it plots all the chords,
                if False, it sees what was deactivated in self.positions and ignores them
            positive_color: color for the positive chords. Default: blue
            negative_color: color for the negative chords. Default: green
            error_color: color for the error bars of the chords. Default: red
        """
        positions = self.positions
        for site in positions.keys():
            pos_obs = positions[site]
            for lc in pos_obs.keys():
                pos_lc = pos_obs[lc]
                if type(pos_lc) != _PositionDict:
                    continue
                if pos_lc['status'] == 'negative':
                    arr = np.array([pos_lc['start_obs']['value'], pos_lc['end_obs']['value']])
                    plt.plot(*arr.T, '--', color=negative_color, linewidth=1.5)
                else:
                    n = 0
                    if pos_lc['immersion']['on'] or all_chords:
                        arr = np.array([pos_lc['immersion']['error']])
                        plt.plot(*arr.T, color=error_color, linewidth=3.0)
                        n+=1
                    if pos_lc['emersion']['on'] or all_chords:
                        arr = np.array([pos_lc['emersion']['error']])
                        plt.plot(*arr.T, color=error_color, linewidth=3.0)
                        n+=1
                    if n == 2:
                        arr = np.array([pos_lc['immersion']['value'], pos_lc['emersion']['value']])
                        plt.plot(*arr.T, color=positive_color, linewidth=2.0)
        plt.axis('equal')

    def get_map_sites(self):
        """Return Dictionary with sites in the format required by plot_occ_map function

        Returns:
            sites (dict): Dictionary with the sites in the format required by plot_occ_map function
        """
        sites = {}
        color = {'positive': 'blue', 'negative': 'red'}
        for o, l in self.__observations:
            sites[o.name] = [o.lon.deg, o.lat.deg, 10, 10, color[self.positions[o.name][l.name]['status']]]
        return sites

    def plot_occ_map(self, **kwargs):
        """Plot occultation map

        Parameters:
            All parameters are parsed directly by prediction.plot_occ_map()
            Please refer to the tutorial
        """
        kwargs['sites'] = kwargs.get('sites', self.get_map_sites())
        if 'offset' not in kwargs and hasattr(self, 'fitted_params'):
            off_ra = self.fitted_params['center_f'][0]*u.km
            off_dec = self.fitted_params['center_g'][0]*u.km
            off_ra = np.arctan2(off_ra, self.dist)
            off_dec = np.arctan2(off_dec, self.dist)
            kwargs['offset'] = [off_ra, off_dec]
        self.predict.plot_occ_map(**kwargs)

    def __str__(self):
        """String representation of the Star class
        """
        out = 'Stellar occultation of star Gaia-DR2 {} by {}.\n\n'.format(self.star.code, self.ephem.name)
        out += 'Geocentric Closest Approach: {:.3f}\n'.format(self.ca)
        out += 'Instant of CA: {}\n'.format(self.tca.iso)
        out += 'Position Angle: {:.2f}\n'.format(self.pa)
        out += 'Geocentric shadow velocity: {:.2f}\n\n\n'.format(self.vel)

        count = {'positive': 0, 'negative': 0, 'visual': 0}
        string = {'positive': '', 'negative': '', 'visual': ''}
        try:
            pos = self.positions
        except:
            pass
        for o,l in self.__observations:
            status = pos[o.name][l.name]['status']
            if count[status] > 0:
                string[status] += '-'*79 + '\n'
            string[status] += o.__str__() + '\n\n'
            string[status] += l.__str__() + ''
            count[status] += 1

        if np.sum([count[i] for i in count.keys()]) == 0:
            out += 'No observations reported'
        else:
            out += '\n'.join(['{} {} observations'.format(count[k],k) for k in string.keys() if count[k] > 0])

        out += '\n\n'

        out += '#'*79 + '\n{:^79s}\n'.format('STAR') + '#'*79 + '\n'
        out += self.star.__str__() + '\n\n'

        coord = self.star.geocentric(self.tca)
        error_star = self.star.error_at(self.tca)
        out += 'GCRS star coordinate at occultation Epoch ({}):\n'.format(self.tca.iso)
        out += 'RA={} +/- {:.4f}, DEC={} +/- {:.4f}\n\n'.format(
            coord.ra.to_string(u.hourangle, sep='hms', precision=5), error_star[0],
            coord.dec.to_string(u.deg, sep='dms', precision=4), error_star[1])
        
        out += '#'*79 + '\n{:^79s}\n'.format('EPHEMERIS') + '#'*79 +'\n'
        out += self.ephem.__str__() + '\n'

        for status in string.keys():
            if count[status] == 0:
                continue
            out += ('#'*79 + '\n{:^79s}\n'.format(status.upper() + ' OBSERVATIONS') + '#'*79)
            out += '\n\n'
            out += string[status]

        if hasattr(self, 'fitted_params'):
            out += '#'*79 + '\n{:^79s}\n'.format('RESULTS') + '#'*79 + '\n\n'
            out += 'Fitted Ellipse:\n'
            out += '\n'.join(['{}: {:.3f} +/- {:.3f}'.format(k, *self.fitted_params[k])
                              for k in self.fitted_params.keys()]) + '\n'
            out += '\nMinimum chi-square: {:.3f}\n'.format(self.chi2_params['chi2_min'])
            out += 'Number of fitted points: {}\n'.format(self.chi2_params['npts'])
            out += 'Number of fitted parameters: {}\n'.format(self.chi2_params['nparam'])
            out += 'Radial dispersion: {:.3f} km\n'.format(self.chi2_params['radial_dispersion'])
            out += 'Mean error: {:.3f} +/- {:.3f} km\n'.format(self.chi2_params['mean_error'][0], self.chi2_params['mean_error'][1])
            
            out += '\n' + self.new_astrometric_position(log=False)

        return out
