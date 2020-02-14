import numpy as np
from astropy.coordinates import SkyCoord, SkyOffsetFrame
import astropy.units as u
from .config import test_attr
from .star import Star

class EphemPlanete():
    """ EphemPlanete simulates ephem_planete and fit_d2_ksi_eta.

    Parameters:
        ephem (str):Input file with hour, minute, RA, DEC, distance  

    """
    def __init__(self, ephem):
        data = np.loadtxt(ephem, unpack=True)
        self.time = data[0]
        self.__reftime = self.time[0]
        self.ephem = SkyCoord(data[1]*u.deg, data[2]*u.deg, data[3]*u.AU)
        self.min_time = self.time.min()
        self.max_time = self.time.max()

    def fit_d2_ksi_eta(self, star):
        """ Fits the on-sky ephemeris position relative to a star
        
        Parameters:
        star (str, SkyCoord):The coordinate of the star in the same frame as the ephemeris.
        """
        if hasattr(self, 'star') and self.star == star:
            return
        if type(star) == str:
            star = SkyCoord(star, unit=(u.hourangle, u.deg))
        self.star = star
        target = self.ephem.transform_to(SkyOffsetFrame(origin=star))  
        da = -target.cartesian.y
        dd = -target.cartesian.z
        dt = (self.time-self.__reftime)/(self.max_time-self.min_time)

        self.ksi = np.polyfit(dt, da.to(u.km).value, 2)
        self.eta = np.polyfit(dt, dd.to(u.km).value, 2)
        
    def get_ksi_eta(self, time, star=None):
        """ Returns the on-sky position of the ephemeris relative to a star.
        
        Parameters:
        time (int, float):Time from which to calculate the position.
        star (str, SkyCoord):The coordinate of the star in the same frame as the ephemeris.
        
        Returns:
        ksi, eta (float): on-sky position of the ephemeris relative to a star
        """
        if star:
            self.fit_d2_ksi_eta(star)
        if time < self.min_time or time > self.max_time:
            raise ValueError('time must be in the interval [{},{}]'.format(self.min_time, self.max_time))
        if hasattr(self, 'ksi') and hasattr(self, 'eta'):
            ksi = np.poly1d(self.ksi)
            eta = np.poly1d(self.eta)
            return ksi((time-self.__reftime)/(self.max_time-self.min_time)), eta((time-self.__reftime)/(self.max_time-self.min_time))
        else:
            raise ValueError('A "star" parameter is missing. Please run fit_d2_ksi_eta first.')

    def __str__(self):
        """ Print the ksi, eta fit values
        """
        out = 'Ephemeris valid from {} until {}\n'.format(self.min_time, self.max_time)
        if hasattr(self, 'ksi') and hasattr(self, 'eta'):
            out += 'star coordinates: {}\n'.format(self.star.to_string('hmsdms'))
            out += 'fit for time=(jd-{})/({}-{})\n\n'.format(self.__reftime, self.max_time,self.min_time)
            out += '                 aksi={}\n'.format(self.ksi[0])
            out += '                 bksi={}\n'.format(self.ksi[1])
            out += '                 cksi={}\n'.format(self.ksi[2])
            out += '                 aeta={}\n'.format(self.eta[0])
            out += '                 beta={}\n'.format(self.eta[1])
            out += '                 ceta={}\n'.format(self.eta[2])
        return out

### Object for ephemeris
class Ephemeris(EphemPlanete):
    '''
    Docstring
    Define ephemeris
    It can be the style of ephem_planete, ephemeris from JPL or with bsp files
    '''
    pass
    #def __init__(self, **kwargs):
    #    if 'ephem' in kwargs:
    #        try:
    #            ephem = test_attr(kwargs['ephem'], str, 'ephem')
    #            self.ephem = EphemPlanete(ephem)
    #        except:
    #            pass
    #    else:
    #        raise InputError('Input values does not correspont to any allowed value')
            
    #def get_position(self, time):
        # returns the position for a given time, it can return ksi, eta
        #time = test_attr(time, Time, 'time')

        #if self.data:
        #    ephem_frame = SkyOffsetFrame(origin=pos)
        #    new_pos = SkyCoord(lon=self.delta.d_lon_coslat, lat=self.delta.d_lat, frame=star_frame)
        #    return new_pos.transform_to(ICRS)
        #return new_pos
    #    return
    
    #def get_topocentric(self, site):
        # return topocentric position given site
    #    return
    
    #def get_ksi_eta(self, time, star=None):
    #    if type(self.ephem) == EphemPlanete:
    #        return self.ephem.get_ksi_eta(time,star)
        # returns the relative position between the ephemeris for a given time and the star
        #star = test_attr(star, Star, 'star')
        #pos = self.get_in(time)
        #target = pos.transform_to(SkyOffsetFrame(origin=star)) 
        #return -target.cartesian.y, -target.cartesian.z
    #    return
    
    #def add_offset(self, da_cosdec, dded):
        # saves an offset for the ephemeris
        #dadc = test_attr(da_cosdec, u.quantity.Quantity, 'd_lon_coslat')
        #dd = test_attr(dded, u.quantity.Quantity, 'dd')
        #self.delta = SphericalCosLatDifferential(dadc, dd, 0.0*u.km)
    #    return
        
    #def __str__(self):
        # return what it is to be printed
    #    return ''
        