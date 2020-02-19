from .config import test_attr
from .star import Star
from .ephem import Ephemeris
from .observer import Observer
import astropy.units as u

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
    if type(ephem) != Ephemeris:
        raise ValueError('ephem must be a Ephemeris object')
    if type(observer) != Observer:
        raise ValueError('observer must be an Observer object')
        
    coord = star.geocentric(time)
    dt = 0.1*u.s
    
    ksio1, etao1 = observer.get_ksi_eta(time=time, star=coord)
    ksie1, etae1 = ephem.get_ksi_eta(time=time.jd, star=coord)
    
    f = ksio1+ksie1
    g = etao1+etae1
    
    ksio2, etao2 = observer.get_ksi_eta(time=time+dt, star=coord)
    ksie2, etae2 = ephem.get_ksi_eta(time=(time+dt).jd, star=coord)
    
    nf = ksio2+ksie2
    ng = etao2+etae2
    
    vf = (nf-f)/0.1
    vg = (ng-g)/0.1
    ## falta calcular velocidade
    return f, g, vf, vg
        
### Object for occultation
class Occultation():
    '''
    Docstring
    Do the reduction of the occultation
    '''
    def __init__(self, star, ephem):
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
    
    def add_observation(self, obs):
        """ Add Observers to the Occultation object.
        
        Parameters:
        obs (Observer):The Observer object to be added.

        """
        if type(star) != Observer:
            raise ValueError('obs must be an Observer object')
        if not hasattr(self, obs):
            self.obs = []
        self.obs.append(obs)
    
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
        # return what it is to be printed
        return ''