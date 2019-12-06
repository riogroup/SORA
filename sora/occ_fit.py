from .config import test_attr
from .star import Star
from .ephem import Ephem
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
        # Run initial parameters
        #self.star = test_attr(star, Star, 'star')
        #self.ephem = test_attr(ephem, Ephem, 'star')
        return
    
    def add_observation(self, obs):
        # add an observation object to the list of observations of this occultation
        #if not hasattr(self, obs):
        #    self.obs = []
        #self.obs.append(test_attr(obs, Observer, 'Occ_obs'))
        return
    
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