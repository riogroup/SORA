from astropy.coordinates import SkyCoord ICRS, SkyOffsetFrame, SphericalCosLatDifferential
from astropy.time import Time
import astropy.units as u
import spiceypy as spice
from astroquery.vizier import Vizier
from .config import test_attr
from .star import Star


### Object for ephemeris
class Ephem():
    '''
    Docstring
    Define ephemeris
    It can be the style of ephem_planete, ephemeris from JPL or with bsp files
    '''
    def __init__(self, **kwargs):
        # Run initial parameters
        #if 'ephem' in kwargs:
        #    self.ephem = test_attr(kwargs['ephem'], SkyCoord, 'ephem')
        #if 'name' in kwargs:
        #    self.name = kwargs['name']
        return
            
    def get_position(self, time):
        # returns the position for a given time, it can return ksi, eta
        #time = test_attr(time, Time, 'time')

        #if self.data:
        #    ephem_frame = SkyOffsetFrame(origin=pos)
        #    new_pos = SkyCoord(lon=self.delta.d_lon_coslat, lat=self.delta.d_lat, frame=star_frame)
        #    return new_pos.transform_to(ICRS)
        #return new_pos
        return
    
    def get_relat_position(self, star, time):
        # returns the relative position between the ephemeris for a given time and the star
        #star = test_attr(star, Star, 'star')
        #pos = self.get_in(time)
        #target = pos.transform_to(SkyOffsetFrame(origin=star)) 
        #return -target.cartesian.y, -target.cartesian.z
        return
    
    def add_offset(self, da_cosdec, dded):
        # saves an offset for the ephemeris
        #dadc = test_attr(da_cosdec, u.quantity.Quantity, 'd_lon_coslat')
        #dd = test_attr(dded, u.quantity.Quantity, 'dd')
        #self.delta = SphericalCosLatDifferential(dadc, dd, 0.0*u.km)
        return
        
    def __str__(self):
        # return what it is to be printed
        return ''
        