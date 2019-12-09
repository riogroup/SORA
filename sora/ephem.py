from astropy.coordinates import SkyCoord, ICRS, SkyOffsetFrame, SphericalCosLatDifferential
from astropy.time import Time
import astropy.units as u
import spiceypy as spice
from astroquery.vizier import Vizier
from .config import test_attr
from .star import Star

class EphemPlanete():
    """ EphemPlanete simulates ephem_planete and fit_d2_ksi_eta.

    Parameters:
        input_file (str):Input file with Time

    Returns:
        catalogue(astropy.Table):An astropy Table with the catalogue informations.   

    """
    def __init__(self, input_file):
        data = np.loadtxt(input_file, unpack=True)
        self.time = data[0]+data[1]/60.0
        self.ephem = SkyCoord(a[2]*u.deg, a[3]*u.deg, a[4]*u.AU)
        self.min_time(self.time.min())
        self.max_time(self.time.max())

    def fit_d2_ksi_eta(self, star):
        if hasattr(self, "star") and self.star == star:
            continue
        self.star = star
        target = ephem.transform_to(SkyOffsetFrame(origin=star))  
        da = -target.cartesian.y
        dd = -target.cartesian.z

        self.ksi = np.polyfit(self.time, da.to(u.km).value, 2)
        self.eta = np.polyfit(self.time, dd.to(u.km).value, 2)
        
        


### Object for ephemeris
class Ephemeris():
    '''
    Docstring
    Define ephemeris
    It can be the style of ephem_planete, ephemeris from JPL or with bsp files
    '''
    def __init__(self, **kwargs):
        if 'ephem' in kwargs:
            try:
                ephem = test_attr(kwargs['ephem'], str, 'ephem')
                self.ephem = EphemPlanete(ephem)
            except:
                pass
        else:
            raise InputError('Input values does not correspont to any allowed value')
            
    def get_position(self, time):
        # returns the position for a given time, it can return ksi, eta
        #time = test_attr(time, Time, 'time')

        #if self.data:
        #    ephem_frame = SkyOffsetFrame(origin=pos)
        #    new_pos = SkyCoord(lon=self.delta.d_lon_coslat, lat=self.delta.d_lat, frame=star_frame)
        #    return new_pos.transform_to(ICRS)
        #return new_pos
        return
    
    def get_topocentric(self, site):
        # return topocentric position given site
        return
    
    def get_ksi_eta(self, star, time):
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
        