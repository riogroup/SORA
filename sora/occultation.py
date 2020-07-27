from .star import Star
from .ephem import EphemPlanete, EphemJPL, EphemKernel
from .observer import Observer
from .lightcurve import LightCurve
from .prediction import occ_params, PredictionTable
from .extra import ChiSquare
from sora.config.decorators import deprecated_alias
import astropy.units as u
from astropy.coordinates import SkyCoord, SkyOffsetFrame
from astropy.time import Time
import numpy as np
import warnings
import matplotlib.pyplot as plt


warnings.simplefilter('always', UserWarning)


def positionv(star, ephem, observer, time):
    """ Calculates the position and velocity of the occultation shadow relative to the observer.

    Parameters:
        star (Star): The coordinate of the star in the same reference frame as the ephemeris.
           It must be a Star object.
        ephem (Ephem): The object ephemeris. It must be an Ephemeris object.
        observer (Observer): The Observer information. It must be an Observer object.
        time (Time): Reference instant to calculate position and velocity.

    Return:
        f, g (float): The orthographic projection of the shadow relative to the observer.
            f is in the x-axis (East-West direction; East positive)
            g is in the y-axis (North-South direction; North positive)
    """
    if type(star) != Star:
        raise ValueError('star must be a Star object')
    if type(ephem) not in [EphemPlanete, EphemJPL, EphemKernel]:
        raise ValueError('ephem must be an Ephemeris object')
    if type(observer) != Observer:
        raise ValueError('observer must be an Observer object')
    time = Time(time)

    coord = star.geocentric(time)
    dt = 0.1*u.s

    if type(ephem) == EphemPlanete:
        ephem.fit_d2_ksi_eta(coord, log=False)
    ksio1, etao1 = observer.get_ksi_eta(time=time, star=coord)
    ksie1, etae1 = ephem.get_ksi_eta(time=time, star=coord)

    f = ksio1-ksie1
    g = etao1-etae1

    ksio2, etao2 = observer.get_ksi_eta(time=time+dt, star=coord)
    ksie2, etae2 = ephem.get_ksi_eta(time=time+dt, star=coord)

    nf = ksio2-ksie2
    ng = etao2-etae2

    vf = (nf-f)/0.1
    vg = (ng-g)/0.1

    return f, g, vf, vg


@deprecated_alias(pos_angle='position_angle', dpos_angle='dposition_angle')  # remove this line for v1.0
def fit_ellipse(*args, equatorial_radius, dequatorial_radius=0, center_f=0, dcenter_f=0, center_g=0,
                dcenter_g=0, oblateness=0, doblateness=0, position_angle=0, dposition_angle=0,
                loop=10000000, number_chi=10000, dchi_min=None, log=False):
    """ Fits an ellipse to given occultation using given parameters

    Parameters:
        required params:
        Each occultation is added as the first arguments directly.
            center_f (int,float): The coordinate in f of the ellipse center. Default=0
            center_g (int,float): The coordinate in g of the ellipse center. Default=0
            equatorial_radius (int,float): The Equatorial radius (semi-major axis) of the ellipse.
            oblateness (int,float): The oblateness of the ellipse. Default=0 (circle)
            position_angle (int,float): The pole position angle of the ellipse in degrees. Default=0
                Zero is in the North direction ('g-positive'). Positive clockwise.

        Parameters interval of fitting. Default values are set to zero.
        Search between (value - dvalue) and (value + dvalue):
            dcenter_f (int,float): Interval for coordinate f of the ellipse center
            dcenter_g (int,float): Interval for coordinate g of the ellipse center
            dequatorial_radius (int,float): Interval for the Equatorial radius (semi-major axis) of the ellipse
            doblateness (int,float): Interval for the oblateness of the ellipse
            dposition_angle (int,float): Interval for the pole position angle of the ellipse in degrees

        loop (int): The number of ellipses to attempt fitting. Default: 10,000,000
        dchi_min (intt,float): If given, it will only save ellipsis which chi square are
            smaller than chi_min + dchi_min.
        number_chi (int): if dchi_min is given, the procedure is repeated until
            number_chi is reached. Default: 10,000
        log (bool): If True, it prints information while fitting. Default: False.

    Returns:
        chisquare: A ChiSquare object with all parameters.

    Examples:
        fit_ellipse(occ1, **kwargs) to fit the ellipse to the chords of occ1 Occultation object
        fit_ellipse(occ1, occ2, **kwargs) to fit the ellipse to the chords of occ1 and occ2 Occultation objects together
    """
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
                        f, g = pos_lc['immersion']['value']
                        err = np.array(pos_lc['immersion']['error'])
                        erro = np.linalg.norm(err[0]-err[1])/2.0
                        values.append([f, g, erro])
                    if pos_lc['emersion']['on']:
                        f, g = pos_lc['emersion']['value']
                        err = np.array(pos_lc['emersion']['error'])
                        erro = np.linalg.norm(err[0]-err[1])/2.0
                        values.append([f, g, erro])

    controle_f0 = Time.now()
    f0_chi = np.array([])
    g0_chi = np.array([])
    a_chi = np.array([])
    obla_chi = np.array([])
    posang_chi = np.array([])
    chi2_best = np.array([])

    while (len(f0_chi) < number_chi):
        chi2 = np.zeros(loop)
        f0 = center_f + dcenter_f*(2*np.random.random(loop) - 1)
        g0 = center_g + dcenter_g*(2*np.random.random(loop) - 1)
        a = equatorial_radius + dequatorial_radius*(2*np.random.random(loop) - 1)
        obla = oblateness + doblateness*(2*np.random.random(loop) - 1)
        obla[obla < 0], obla[obla > 1] = 0, 1
        phi_deg = position_angle + dposition_angle*(2*np.random.random(loop) - 1)
        controle_f1 = Time.now()

        for fi, gi, si in values:
            b = a - a*obla
            phi = phi_deg*(np.pi/180.0)
            dfi = fi-f0
            dgi = gi-g0
            r = np.sqrt(dfi**2 + dgi**2)
            theta = np.arctan2(dgi, dfi)
            ang = theta+phi
            r_model = (a*b)/np.sqrt((a*np.sin(ang))**2 + (b*np.cos(ang))**2)
            f_model = f0 + r_model*np.cos(theta)
            g_model = g0 + r_model*np.sin(theta)
            chi2 += ((fi - f_model)**2 + (gi - g_model)**2)/(si**2)

        controle_f2 = Time.now()
        if dchi_min is not None:
            region = np.where(chi2 < chi2.min() + dchi_min)[0]
        else:
            region = np.arange(len(chi2))
        chi2_best = np.append(chi2_best, chi2[region])
        if log:
            print('Elapsed time: {:.3f} seconds.'.format((controle_f2 - controle_f1).sec))
            print(len(chi2[region]), len(chi2_best))
        f0_chi = np.append(f0_chi, f0[region])
        g0_chi = np.append(g0_chi, g0[region])
        a_chi = np.append(a_chi, a[region])
        obla_chi = np.append(obla_chi, obla[region])
        posang_chi = np.append(posang_chi, phi_deg[region])

    chisquare = ChiSquare(chi2_best, len(values), center_f=f0_chi, center_g=g0_chi, equatorial_radius=a_chi,
                          oblateness=obla_chi, position_angle=posang_chi)
    controle_f4 = Time.now()
    if log:
        print('Total elapsed time: {:.3f} seconds.'.format((controle_f4 - controle_f0).sec))

    onesigma = chisquare.get_nsigma(sigma=1)
    a = onesigma['equatorial_radius'][0]
    f0 = onesigma['center_f'][0]
    g0 = onesigma['center_g'][0]
    obla = onesigma['oblateness'][0]
    phi_deg = onesigma['position_angle'][0]
    radial_dispersion = np.array([])
    error_bar = np.array([])

    for fi, gi, si in values:
        b = a - a*obla
        phi = phi_deg*(np.pi/180.0)
        dfi = fi-f0
        dgi = gi-g0
        r = np.sqrt(dfi**2 + dgi**2)
        theta = np.arctan2(dgi, dfi)
        ang = theta+phi
        r_model = (a*b)/np.sqrt((a*np.sin(ang))**2 + (b*np.cos(ang))**2)
        f_model = f0 + r_model*np.cos(theta)
        g_model = g0 + r_model*np.sin(theta)
        radial_dispersion = np.append(radial_dispersion, r - r_model)
        error_bar = np.append(error_bar, si)

    for occ in args:
        if type(occ) == Occultation:
            occ.fitted_params = {i: onesigma[i] for i in ['equatorial_radius', 'center_f', 'center_g',
                                                          'oblateness', 'position_angle']}
            occ.chi2_params = {'radial_dispersion': [radial_dispersion.mean(), radial_dispersion.std(ddof=1)]}
            occ.chi2_params['mean_error'] = [error_bar.mean(), error_bar.std()]
            occ.chi2_params['chi2_min'] = chisquare.get_nsigma()['chi2_min']
            occ.chi2_params['nparam'] = chisquare.nparam
            occ.chi2_params['npts'] = chisquare.npts
    return chisquare


class _PositionDict(dict):
    """ This is a modified Dictionary object to allow switching on/off of data points.
        It also avoids user to change data.
    """
    def __setitem__(self, key, value):
        """ Redefines how to set a value to a key in the dictionary.
            It only sets a value if the key starts with '_occ_'. Otherwise, it only allows for the user to provide
            'on' or 'off' which is passed only to change the 'on' keyword.
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
            elif key == 'on':
                super().__setitem__('on', status[value])
        else:
            if value not in status.keys():
                raise ValueError("Value must be 'on' or 'off' only.")
            for key1 in self.keys():
                if type(self[key1]) == _PositionDict:
                    if key in self[key1].keys():
                        n = 1
                        self[key1][key] = value
        if n == 0:
            raise KeyError('Key "{}" does not exist'.format(key))

    def __str__(self):
        out = '\n' + '\n'.join(['{}: {}'.format(key, self[key]) for key in self.keys()])
        return out.replace('\n', '\n  ')


class Occultation():
    """ Does the reduction of the occultation
    """
    def __init__(self, star, ephem, time):
        """ Instantiates the Occultation object.

        Parameters:
            star (Star): The coordinate of the star in the same reference frame as the ephemeris.
                It must be a Star object.
            ephem (Ephem): object ephemeris. It must be an Ephemeris object.
            time (str, Time): Reference time of the occultation.
                Time does not need to be exact, but needs to be within approximately 10 minutes
                of the occultation closest approach to calculate occultation parameters.
        """
        if type(star) != Star:
            raise ValueError('star must be a Star object')
        if type(ephem) not in [EphemPlanete, EphemKernel, EphemJPL]:
            raise ValueError('ephem must be a Ephemeris object')
        self.star = star
        self.ephem = ephem

        tca, ca, pa, vel, dist = occ_params(star, ephem, time)
        self.ca = ca   # Closest Approach distance
        self.pa = pa   # Position Angle at CA
        self.vel = vel  # Shadow velocity at CA
        self.dist = dist  # object distance at CA
        self.tca = tca   # Instant of CA
        self.star_diam = self.star.apparent_diameter(self.dist, log=False)

        meta = {
            'name': self.ephem.name, 'radius': self.ephem.radius.to(u.km).value,
            'error_ra': self.ephem.error_ra.to(u.mas).value, 'error_dec': self.ephem.error_dec.to(u.mas).value}
        self.predict = PredictionTable(
            time=[tca], coord_star=[self.star.geocentric(tca)],
            coord_obj=[self.ephem.get_position(tca)], ca=[ca.value], pa=[pa.value], vel=[vel.value],
            dist=[dist.value], mag=[self.star.mag['G']], source=[self.star.code], meta=meta)

        self.__observations = []
        self._position = _PositionDict()

    def add_observation(self, obs, lightcurve):
        """ Adds observations to the Occultation object.

        Parameters:
            obs (Observer): The Observer object to be added.
            lightcurve (LightCurve): The LightCurve object to be added
        """
        if type(obs) != Observer:
            raise ValueError('obs must be an Observer object')
        if type(lightcurve) != LightCurve:
            raise ValueError('lightcurve must be a LightCurve object')
        for o, l in self.__observations:
            if l.name == lightcurve.name:
                raise ValueError('{} LightCurve already associated to {} Observer'.format(lightcurve.name, o.name))
        self.__observations.append((obs, lightcurve))
        lightcurve.set_vel(np.absolute(self.vel))
        lightcurve.set_dist(float(self.dist.AU))
        lightcurve.set_star_diam(float(self.star_diam.km))
        try:
            lightcurve.calc_magnitude_drop(mag_star=self.star.mag['G'], mag_obj=self.ephem.apparent_magnitude(self.tca))
        except:
            lightcurve.bottom_flux = 0.0
            warnings.warn('Magnitude drop was not calculated. Using bottom flux as 0.0 instead.')

    def remove_observation(self, key, key_lc=None):
        """ Removes an observation from the Occultation object.

        Parameters:
            key (str): The name given to Observer or LightCurve to remove from the list.
            keylc (str): In the case where repeated names are present for different observations,
                keylc must be given for the name of the LightCurve and key will be used for the name of the Observer.
        """
        rm_list = np.array([])
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
        ko = np.array(ko)
        kl = np.array(kl)
        if not same_key:
            k = ko[np.where(ko == kl)[0]]
            rm_list = np.hstack((rm_list, k))
        else:
            rm_list = np.hstack((rm_list, np.array(kl)))
            if len(kl) > 0 and len(ko) > 0 and kl[0] not in ko:
                raise ValueError("Observation could not univocally be identified, "
                                 "please give parameters for Observer and LightCurve")
            rm_list = np.hstack((rm_list, np.array(ko)))
        rm_list = np.unique(rm_list)
        if len(rm_list) == 0:
            raise ValueError('No observer "{}" and/or lightcurve "{}" was found'.format(key, key_lc))
        list = np.arange(len(self.__observations)).tolist()
        for i in rm_list:
            list.remove(i)
        self.__observations = [self.__observations[item] for item in list]

    def observations(self):
        """ Print all the observations added to the Occultation object
            Pair (Observer, LightCurve)
        """
        for o, l in self.__observations:
            print('Observer= {}, LC: {}'.format(o.name, l.name))

    def fit_ellipse(self, **kwargs):
        """ Fits an ellipse to the chords of this occultation

        Parameters:
        Required parameters:
            Each occultation is added as the first arguments directly.
            center_f (int,float): The coordinate in f of the ellipse center.
            center_g (int,float): The coordinate in g of the ellipse center.
            equatorial_radius (int,float): The Equatorial radius (semi-major axis) of the ellipse.
            oblateness (int,float): The oblateness of the ellipse.
            position_angle (int,float): The pole position angle of the ellipse in degrees. Default=0
                Zero is in the North direction ('g-positive'). Positive clockwise.

        Parameters interval of fitting. Default values are set to zero.
            Search between (value - dvalue) and (value + dvalue):
            dcenter_f (int,float): Interval for coordinate f of the ellipse center
            dcenter_g (int,float): Interval for coordinate g of the ellipse center
            dequatorial_radius (int,float): Interval for the Equatorial radius (semi-major axis) of the ellipse
            doblateness (int,float): Interval for the oblateness of the ellipse
            dposition_angle (int,float): Interval for the pole position angle of the ellipse in degrees

            loop (int): The number of ellipses to attempt fitting. Default: 10,000,000
            dchi_min (int,float): If given, it will only save ellipsis which chi square are
                smaller than chi_min + dchi_min.
            number_chi (int): if dchi_min is given, the procedure is repeated until
                number_chi is reached. Default: 10,000
            log (bool): If True, it prints information while fitting. Default: False.

        Returns:
            chisquare: A ChiSquare object with all parameters.
        """
        chisquare = fit_ellipse(self, **kwargs)
        return chisquare

    @property
    def positions(self):
        """ Calculates the position and velocity for all chords.
            Saves it into an _PositionDict object.
        """
        position = self._position
        if len(self.__observations) == 0:
            raise ValueError('There is no observation defined for this occultation')

        pair = []
        for o, l in self.__observations:
            pair.append((o.name, l.name))

            coord = [o.lon, o.lat, o.height]
            if o.name not in position.keys():
                position['_occ_'+o.name] = _PositionDict(lon=o.lon, lat=o.lat, height=o.height)
                position[o.name]['_occ_lon'] = o.lon
                position[o.name]['_occ_lat'] = o.lat
                position[o.name]['_occ_height'] = o.height
            pos_obs = position[o.name]
            coord2 = [pos_obs['lon'], pos_obs['lat'], pos_obs['height']]
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
                    f1, g1, vf1, vg1 = positionv(self.star, self.ephem, o, l.immersion)
                    obs_im['_occ_time'] = l.immersion
                    obs_im['_occ_value'] = (round(f1, 3), round(g1, 3))
                    obs_im['_occ_vel'] = (round(vf1, 3), round(vg1, 3))
                if not do_err and 'time_err' in obs_im.keys() and obs_im['time_err'] == l.immersion_err:
                    pass
                else:
                    fe1, ge1 = positionv(self.star, self.ephem, o, l.immersion-l.immersion_err*u.s)[0:2]
                    fe2, ge2 = positionv(self.star, self.ephem, o, l.immersion+l.immersion_err*u.s)[0:2]
                    obs_im['_occ_time_err'] = l.immersion_err
                    obs_im['_occ_error'] = ((round(fe1, 3), round(ge1, 3)), (round(fe2, 3), round(ge2, 3)))

            if hasattr(l, 'emersion'):
                if 'emersion' not in pos_lc.keys():
                    pos_lc['_occ_emersion'] = _PositionDict(on=True)
                obs_em = pos_lc['emersion']
                do_err = False
                if samecoord and 'time' in obs_em.keys() and obs_em['time'] == l.emersion:
                    pass
                else:
                    do_err = True
                    f1, g1, vf1, vg1 = positionv(self.star, self.ephem, o, l.emersion)
                    obs_em['_occ_time'] = l.emersion
                    obs_em['_occ_value'] = (round(f1, 3), round(g1, 3))
                    obs_em['_occ_vel'] = (round(vf1, 3), round(vg1, 3))
                if not do_err and 'time_err' in obs_em.keys() and obs_em['time_err'] == l.emersion_err:
                    pass
                else:
                    fe1, ge1 = positionv(self.star, self.ephem, o, l.emersion-l.emersion_err*u.s)[0:2]
                    fe2, ge2 = positionv(self.star, self.ephem, o, l.emersion+l.emersion_err*u.s)[0:2]
                    obs_em['_occ_time_err'] = l.emersion_err
                    obs_em['_occ_error'] = ((round(fe1, 3), round(ge1, 3)), (round(fe2, 3), round(ge2, 3)))

            if pos_lc['status'] == 'negative':
                if 'start_obs' not in pos_lc.keys():
                    pos_lc['_occ_start_obs'] = _PositionDict(on=True)
                obs_start = pos_lc['start_obs']
                if samecoord and 'time' in obs_start.keys() and obs_start['time'] == l.initial_time:
                    pass
                else:
                    f, g, vf, vg = positionv(self.star, self.ephem, o, l.initial_time)
                    obs_start['_occ_time'] = l.initial_time
                    obs_start['_occ_value'] = (round(f, 3), round(g, 3))
                    obs_start['_occ_vel'] = (round(vf, 3), round(vg, 3))
                if 'end_obs' not in pos_lc.keys():
                    pos_lc['_occ_end_obs'] = _PositionDict(on=True)
                obs_end = pos_lc['end_obs']
                if samecoord and 'time' in obs_end.keys() and obs_end['time'] == l.end_time:
                    pass
                else:
                    f, g, vf, vg = positionv(self.star, self.ephem, o, l.end_time)
                    obs_end['_occ_time'] = l.end_time
                    obs_end['_occ_value'] = (round(f, 3), round(g, 3))
                    obs_end['_occ_vel'] = (round(vf, 3), round(vg, 3))

        for key in list(position):
            n = 0
            for key_lc in list(position[key]):
                if type(position[key][key_lc]) != _PositionDict:
                    continue
                if (key, key_lc) not in pair:
                    del position[key][key_lc]
                else:
                    n += 1
            if n == 0:
                del position[key]

        return self._position

    @positions.setter
    def positions(self, value):
        """ If the users tries to set a value to position, it must be 'on' or 'off',
            and it will be assigned to all chords.
        """
        if value not in ['on', 'off']:
            raise ValueError("Value must be 'on' or 'off' only.")
        pos = self.positions
        for key in pos.keys():
            pos[key] = value

    def check_velocities(self):
        """ Prints the current velocity used by the LightCurves and its Radial velocity.
        """
        if hasattr(self, 'fitted_params'):
            center = np.array([self.fitted_params['center_f'][0], self.fitted_params['center_g'][0]])
        else:
            center = np.array([0, 0])
        positions = self.positions
        for ob, lc in self.__observations:
            vals = positions[ob.name][lc.name]
            if all([i not in vals.keys() for i in ['immersion', 'emersion']]):
                continue
            print('{} - Velocity used: {:.3f}'.format(lc.name, lc.vel))
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
        """ Calculates the new astrometric position for the object given fitted parameters

        Parameters:
            time (str,Time): Reference time to calculate the position.
                If not given, it uses the instant of the occultation Closest Approach.
            offset (list): Offset to apply to the position. If not given, uses the parameters from the fitted ellipse.
                Must be a list of 3 values being [X, Y, 'unit']
                    'unit' must be Angular or Distance unit.
                    - If Distance units for X and Y
                        Ex: [100, -200, 'km'], [0.001, 0.002, 'AU']
                    - If Angular units fox X [d*a*cos(dec)] and Y [d*dec]
                        Ex: [30.6, 20, 'mas'], or [-15, 2, 'arcsec']

            error (list): Error bar of the given offset. If not given, it uses the 1-sigma value of the fitted ellipse.
                Error must be a list of 3 values being [dX, dY, 'unit'], similar to offset.
                It does not need to be in the same unit as offset.
            log (bool): If true, it Prints text, else it Returns text.
        """
        if time is not None:
            time = Time(time)
        else:
            time = self.tca

        tca_diff = np.absolute(time-self.tca)
        if tca_diff > 1*u.day:
            warnings.warn('The difference between the given time and the closest approach instant is {:.1f} days. '
                          'This position could not have a physical meaning.'.format(tca_diff.jd))

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
        error_ra = np.sqrt(error_star[0]**2 + e_off_ra**2)
        error_dec = np.sqrt(error_star[1]**2 + e_off_dec**2)

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

    def plot_chords(self, all_chords=True, positive_color='blue', negative_color='green', error_color='red',
                    ax=None, lw=2, axis_labels=True):
        """ Plots the chords of the occultation

        Parameters:
            all_chords (bool): if True, it plots all the chords,
                if False, it sees what was deactivated in self.positions and ignores them
            positive_color (str): color for the positive chords. Default: blue
            negative_color (str): color for the negative chords. Default: green
            error_color (str): color for the error bars of the chords. Default: red
            ax (maptlotlib.Axes): Axis where to plot chords. Default: Use matplotlib pool.
            lw (int, float): linewidth of the chords. Default: 2
            axis_labels (bool): If True it prints the labels of the axis of the image.
        """
        ax = ax or plt.gca()
        positions = self.positions
        for site in positions.keys():
            pos_obs = positions[site]
            for lc in pos_obs.keys():
                pos_lc = pos_obs[lc]
                if type(pos_lc) != _PositionDict:
                    continue
                if pos_lc['status'] == 'negative':
                    arr = np.array([pos_lc['start_obs']['value'], pos_lc['end_obs']['value']])
                    ax.plot(*arr.T, '--', color=negative_color, linewidth=lw)
                else:
                    n = 0
                    if pos_lc['immersion']['on'] or all_chords:
                        arr = np.array([pos_lc['immersion']['error']])
                        ax.plot(*arr.T, color=error_color, linewidth=lw*3/2)
                        n += 1
                    if pos_lc['emersion']['on'] or all_chords:
                        arr = np.array([pos_lc['emersion']['error']])
                        ax.plot(*arr.T, color=error_color, linewidth=lw*3/2)
                        n += 1
                    if n == 2:
                        arr = np.array([pos_lc['immersion']['value'], pos_lc['emersion']['value']])
                        ax.plot(*arr.T, color=positive_color, linewidth=lw)
        ax.invert_xaxis()
        if axis_labels:
            ax.set_xlabel('f (km)')
            ax.set_ylabel('g (km)')
        ax.axis('equal')

    def get_map_sites(self):
        """ Returns Dictionary with sites in the format required by plot_occ_map function

        Returns:
            sites (dict): Dictionary with the sites in the format required by plot_occ_map function
        """
        sites = {}
        color = {'positive': 'blue', 'negative': 'red'}
        for o, l in self.__observations:
            sites[o.name] = [o.lon.deg, o.lat.deg, 10, 10, color[self.positions[o.name][l.name]['status']]]
        return sites

    def plot_occ_map(self, **kwargs):
        """ Plots the occultation map

        Parameters:
            radius: The radius of the shadow. If not given it uses the equatorial radius
                from the ellipse fit, else it uses the radius obtained from ephem upon instantiating.
            nameimg (str): Change the name of the image saved.
            path (str): Path to a directory where to save map.
            resolution (int): Cartopy feature resolution. "1" means a resolution of "10m",
                "2" a resolution of "50m" and "3" a resolution of "100m". Default = 2
            states (bool): True to plot the state division of the countries. The states of
                some countries will only be shown depending on the resolution.
            zoom (int, float): Zooms in or out of the map.
            centermap_geo (list): Center the map for a given coordinates in longitude and latitude.
                It must be a list with two numbers. Default=None.
            centermap_delta (list): Displace the center of the map given displacement
                in X and Y, in km. It must be a list with two numbers. Default=None.
            centerproj (list): Rotates the Earth to show occultation with the center
                projected at a given longitude and latitude. It must be a list with two numbers
            labels (bool): Plots text above and below the map with the occultation parameters.
                Default=True.
            meridians (int): Plots lines representing the meridians for given interval. Default=30 deg
            parallels (int): Plots lines representing the parallels for given interval. Default=30 deg
            sites (dict): Plots site positions in map. It must be a python dictionary where the key is
                the name of the site, and the value is a list with longitude, latitude, delta_x,
                delta_y and color. delta_x and delta_y are displacement, in km, from the point
                of the site in the map and the name. color is the color of the point.
                If not given, it calculates from observations added to Occultation
            site_name (bool): If True, it prints the name of the sites given, else it plots only the points
            countries (dict): Plots the names of countries. It must be a python dictionary where the key
                is the name of the country and the value is a list with longitude and latitude
                of the lower left part of the text.
            offset (list): applies an offset to the ephemeris, calculating new CA and instant of CA.
                It is a pair of delta_RA*cosDEC and delta_DEC.
                If not given it uses the center from ellipse fitted.
            mapstyle (int): Define the color style of the map. 1 is the default black and white scale.
                "2" is a colored map.
            error (int,float): Ephemeris error in mas. It plots a dashed line representing radius + error.
            ercolor (str): Changes the color of the lines of the error bar.
            ring (int,float): It plots a dashed line representing the location of a ring.
                It is given in km, from the center.
            rncolor (str): Changes the color of ring lines.
            atm (int,float): plots a dashed line representing the location of an atmosphere.
                It is given in km, from the center.
            atcolor (str): Changes the color of atm lines.
            chord_delta (list): list with distances from center to plot chords
            chord_geo (2d-list): list with pairs of coordinates to plot chords
            chcolor (str): color of the line of the chords. Default: grey
            heights (list): plots a circular dashed line showing the locations where the observer
                would observe the occultation at a given height above the horizons.
                This must be a list.
            hcolor (str): Changes the color of the height lines.
            mapsize (list): The size of figure, in cm. It must be a list with two values.
                Default = [46.0, 38.0].
            cpoints (int,float): Interval for the small points marking the center of shadow,
                in seconds. Default=60.
            ptcolor (str): Change the color of the center points.
            alpha (float): The transparency of the night shade, where 0.0 is full transparency
                and 1.0 is full black. Default = 0.2.
            fmt (str): The format to save the image. It is parsed directly by matplotlib.pyplot.
                Default = 'png'
            dpi (int): "Dots per inch". It defines the quality of the image. Default = 100.
            lncolor (str): Changes the color of the line that represents the limits of the shadow over Earth.
            outcolor (str): Changes the color of the lines that represents the limits of the shadow outside Earth
            nscale (int,float): Arbitrary scale for the size of the name of the site.
            cscale (int,float): Arbitrary scale for the name of the country.
            sscale (int,float): Arbitrary scale for the size of point of the site.
            pscale (int,float): Arbitrary scale for the size of the points that represent the center of the shadow
            arrow (bool): If true, it plots the arrow with the occultation direction.

            Comment: Only one of centermap_geo and centermap_delta can be given
        """
        if 'radius' not in kwargs and hasattr(self, 'fitted_params'):
            r_equa = self.fitted_params['equatorial_radius'][0]
            obla = self.fitted_params['oblateness'][0]
            pos_ang = self.fitted_params['position_angle'][0]
            theta = np.linspace(-np.pi, np.pi, 1800)
            map_pa = self.predict['P/A'][0]
            circle_x = r_equa*np.cos(theta)
            circle_y = r_equa*(1.0-obla)*np.sin(theta)
            ellipse_y = -circle_x*np.sin((pos_ang-map_pa)*u.deg) + circle_y*np.cos((pos_ang-map_pa)*u.deg)
            kwargs['radius'] = ellipse_y.max()
            print('Projected shadow radius = {:.1f} km'.format(kwargs['radius']))
        kwargs['sites'] = kwargs.get('sites', self.get_map_sites())
        if 'offset' not in kwargs and hasattr(self, 'fitted_params'):
            off_ra = self.fitted_params['center_f'][0]*u.km
            off_dec = self.fitted_params['center_g'][0]*u.km
            off_ra = np.arctan2(off_ra, self.dist)
            off_dec = np.arctan2(off_dec, self.dist)
            kwargs['offset'] = [off_ra, off_dec]
        self.predict.plot_occ_map(**kwargs)

    def to_log(self, namefile=None):
        """ Saves the occultation log to a file

        Parameters:
            namefile (str): Filename to save the log
        """
        if namefile is None:
            namefile = 'occ_{}_{}.log'.format(self.ephem.name, self.tca.isot[:16])
        f = open(namefile, 'w')
        f.write(self.__str__())
        f.close()

    def to_file(self):
        """ Saves the occultation data to a file

        Three files are saved containing the positions and velocities for the observations.
        They are for the positive, negative and error bars positions.

        The format of the files are: positions in f and g, velocities in f and g, the Julian Date of the observation,
        light curve name of the corresponding position.
        """
        positions = self.positions
        pos = []
        neg = []
        err = []
        for ob, lc in self.__observations:
            status = positions[ob.name][lc.name]['status']
            l_name = lc.name.replace(' ', '_')
            if status == 'positive':
                im = positions[ob.name][lc.name]['immersion']
                f, g = im['value']
                err1, err2 = im['error']
                f1, g1 = err1
                f2, g2 = err2
                vf, vg = im['vel']
                pos.append([f, g, vf, vg, im['time'].jd, l_name+'_immersion'])
                err.append([f1, g1, vf, vg, (im['time']-im['time_err']*u.s).jd, l_name+'_immersion_err-'])
                err.append([f2, g2, vf, vg, (im['time']+im['time_err']*u.s).jd, l_name+'_immersion_err+'])

                em = positions[ob.name][lc.name]['emersion']
                f, g = em['value']
                err1, err2 = em['error']
                f1, g1 = err1
                f2, g2 = err2
                vf, vg = em['vel']
                pos.append([f, g, vf, vg, em['time'].jd, l_name+'_emersion'])
                err.append([f1, g1, vf, vg, (em['time']-em['time_err']*u.s).jd, l_name+'_emersion_err-'])
                err.append([f2, g2, vf, vg, (em['time']+em['time_err']*u.s).jd, l_name+'_emersion_err+'])
            if status == 'negative':
                ini = positions[ob.name][lc.name]['start_obs']
                f, g = ini['value']
                vf, vg = ini['vel']
                neg.append([f, g, vf, vg, ini['time'].jd, l_name+'_start'])

                end = positions[ob.name][lc.name]['end_obs']
                f, g = end['value']
                vf, vg = end['vel']
                neg.append([f, g, vf, vg, end['time'].jd, l_name+'_end'])
        if len(pos) > 0:
            f = open('occ_{}_pos.txt'.format(self.ephem.name), 'w')
            for line in pos:
                f.write('{:9.3f} {:9.3f} {:-6.2f} {:-6.2f} {:16.8f} {}\n'.format(*line))
            f.close()
            f = open('occ_{}_err.txt'.format(self.ephem.name), 'w')
            for line in err:
                f.write('{:9.3f} {:9.3f} {:-6.2f} {:-6.2f} {:16.8f} {}\n'.format(*line))
            f.close()
        if len(neg) > 0:
            f = open('occ_{}_neg.txt'.format(self.ephem.name), 'w')
            for line in neg:
                f.write('{:9.3f} {:9.3f} {:-6.2f} {:-6.2f} {:16.8f} {}\n'.format(*line))
            f.close()

    def __str__(self):
        """ String representation of the Occultation class
        """
        out = ('Stellar occultation of star Gaia-DR2 {} by {}.\n\n'
               'Geocentric Closest Approach: {:.3f}\n'
               'Instant of CA: {}\n'
               'Position Angle: {:.2f}\n'
               'Geocentric shadow velocity: {:.2f}\n'
               'Sun-Geocenter-Target angle:  {:.2f} deg\n'
               'Moon-Geocenter-Target angle: {:.2f} deg\n\n\n'.format(
                   self.star.code, self.ephem.name, self.ca, self.tca.iso,
                   self.pa, self.vel, self.predict['S-G-T'].data[0],
                   self.predict['M-G-T'].data[0])
               )

        count = {'positive': 0, 'negative': 0, 'visual': 0}
        string = {'positive': '', 'negative': '', 'visual': ''}
        if len(self.__observations) > 0:
            pos = self.positions
            for ob, lc in self.__observations:
                status = pos[ob.name][lc.name]['status']
                if count[status] > 0:
                    string[status] += '-'*79 + '\n'
                string[status] += ob.__str__() + '\n'
                ephem_altaz = ob.altaz(lc.time_mean, self.ephem.get_position(lc.time_mean))
                string[status] += 'Target altitude: {:.1f} deg\n'.format(ephem_altaz[0])
                string[status] += 'Target azimuth:  {:.1f} deg\n\n'.format(ephem_altaz[1])
                string[status] += lc.__str__() + ''
                count[status] += 1

        if np.sum([count[i] for i in count.keys()]) == 0:
            out += 'No observations reported'
        else:
            out += '\n'.join(['{} {} observations'.format(count[k], k) for k in string.keys() if count[k] > 0])

        out += '\n\n'

        out += '#'*79 + '\n{:^79s}\n'.format('STAR') + '#'*79 + '\n'
        out += self.star.__str__() + '\n\n'

        coord = self.star.geocentric(self.tca)
        error_star = self.star.error_at(self.tca)
        out += 'Geocentric star coordinate at occultation Epoch ({}):\n'.format(self.tca.iso)
        out += 'RA={} +/- {:.4f}, DEC={} +/- {:.4f}\n\n'.format(
            coord.ra.to_string(u.hourangle, sep='hms', precision=5), error_star[0],
            coord.dec.to_string(u.deg, sep='dms', precision=4), error_star[1])

        out += '#'*79 + '\n{:^79s}\n'.format('EPHEMERIS') + '#'*79 + '\n'
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
            polar_radius = self.fitted_params['equatorial_radius'][0]*(1.0-self.fitted_params['oblateness'][0])
            equivalent_radius = np.sqrt(self.fitted_params['equatorial_radius'][0]*polar_radius)
            out += 'polar_radius: {:.3f} km \n'.format(polar_radius)
            out += 'equivalent_radius: {:.3f} km \n'.format(equivalent_radius)
            if self.ephem.H is not np.nan:
                H_sun = -26.74
                geometric_albedo = (10**(0.4*(H_sun - self.ephem.H))) * ((u.au.to('km')/equivalent_radius)**2)
                out += 'geometric albedo (V): {:.3f} ({:.1%}) \n'.format(geometric_albedo, geometric_albedo)
            else:
                out += 'geometric albedo (V): not calculated, absolute magnitude (H) is unknown \n'
            out += '\nMinimum chi-square: {:.3f}\n'.format(self.chi2_params['chi2_min'])
            out += 'Number of fitted points: {}\n'.format(self.chi2_params['npts'])
            out += 'Number of fitted parameters: {}\n'.format(self.chi2_params['nparam'])
            out += 'Minimum chi-square per degree of freedom: {:.3f}\n'.format(
                self.chi2_params['chi2_min']/(self.chi2_params['npts'] - self.chi2_params['nparam']))
            out += 'Radial dispersion: {:.3f} +/- {:.3f} km\n'.format(
                self.chi2_params['radial_dispersion'][0], self.chi2_params['radial_dispersion'][1])
            out += 'Mean error: {:.3f} +/- {:.3f} km\n'.format(
                self.chi2_params['mean_error'][0], self.chi2_params['mean_error'][1])

            out += '\n' + self.new_astrometric_position(log=False)

        return out
