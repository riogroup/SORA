import .config as cnfg
from .star import Star
from .ephem import Ephemeris
from .observer import Observer

def positionv(ephem,star,observer):
    # calculates the position and velocity given the ephem, star and observer
    return

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