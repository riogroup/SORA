from .config import test_attr
from .star import Star
from .ephem import Ephemeris, EphemPlanete
from .observer import Observer
import astropy.units as u
from astropy.time import Time
import numpy as np
from astropy.coordinates import SphericalCosLatDifferential, SkyCoord, SkyOffsetFrame

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
    if type(ephem) != Ephemeris:
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
    
    pa = np.arctan2(-ksi[min],-eta[min]).to(u.deg)
    
    dksi = ksi[min+1]-ksi[min]
    deta = eta[min+1]-eta[min]
    vel = np.sqrt(dksi**2 + deta**2)/delta_t
    vel = -vel*np.sign(dksi)*(u.km/u.s)
    
    return tt[min], ca, pa, vel, dist
    
        
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
        if type(ephem) != Ephemeris:
            raise ValueError('ephem must be a Ephemeris object')
        self.star = star
        self.ephem = ephem
        
        tt, ca, pa, vel, dist = occ_params(star,ephem, time)
        self.ca = ca   # Closest Approach distance
        self.pa = pa   # Position Angle at CA
        self.vel = vel  # Shadow velocity at CA
        self.dist = dist  # object distance at CA
        self.tca = tt   # Instant of CA
        
        self.obs_positive = []
        self.obs_negative = []
        self.obs_visual = []
        self.obs_undefined = []
    
    def add_observation(self, obs, status='undefined'):
        """ Add Observers to the Occultation object.
        
        Parameters:
        obs (Observer):The Observer object to be added.
        status (string): it can be "positive", "negative", "visual" or "undefined"

        """
        if type(obs) != Observer:
            raise ValueError('obs must be an Observer object')
        if status not in ['positive', 'negative', 'visual', 'undefined']:
            raise ValueError('status must be one the following valeus: "positive", "negative", "visual" or "undefined"')
        if obs in self.obs_positive:
            raise ValueError('{} observation already defined as positive'.format(obs.name))
        elif obs in self.obs_negative:
            raise ValueError('{} observation already defined as negative'.format(obs.name))
        elif obs in self.obs_visual:
            raise ValueError('{} observation already defined as visual'.format(obs.name))
        elif obs in self.obs_undefined:
            raise ValueError('{} observation already defined as undefined'.format(obs.name))
        if status == 'positive':
            self.obs_positive.append(obs)
        if status == 'negative':
            self.obs_negative.append(obs)
        if status == 'visual':
            self.obs_visual.append(obs)
        if status == 'undefined':
            self.obs_undefined.append(obs)
    
    def fit_ellipse(self):
        # fit ellipse to the points
        return
    
    def fit_to_shape(self):
        # fit points to a 3D shape model
        return
    
    def plot_chords(self):
        # plot chords of the occultation
        return
    
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