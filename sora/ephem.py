import numpy as np
from astropy.coordinates import SkyCoord, SkyOffsetFrame
from astropy.time import Time
import astropy.units as u
from .config import test_attr
from .star import Star

class EphemPlanete():
    """ EphemPlanete simulates ephem_planete and fit_d2_ksi_eta.

    Parameters:
        ephem (str):Input file with hour, minute, RA, DEC, distance  

    """
    def __init__(self, name, ephem):
        data = np.loadtxt(ephem, unpack=True)
        self.name = name
        self.time = Time(data[0], format='jd')
        self.ephem = SkyCoord(data[1]*u.deg, data[2]*u.deg, data[3]*u.AU)
        self.__reftime = self.time[0]
        self.min_time = Time(data[0].min(), format='jd')
        self.max_time = Time(data[0].max(), format='jd')

    def fit_d2_ksi_eta(self, star):
        """ Fits the on-sky ephemeris position relative to a star
        
        Parameters:
        star (str, SkyCoord):The coordinate of the star in the same frame as the ephemeris.
        """
        if type(star) == str:
            star = SkyCoord(star, unit=(u.hourangle, u.deg))
        if hasattr(self, 'star') and self.star.to_string('hmsdms', precision=5) == star.to_string('hmsdms', precision=5):
            return
        self.star = star
        target = self.ephem.transform_to(SkyOffsetFrame(origin=star))  
        da = -target.cartesian.y
        dd = -target.cartesian.z
        dt = (self.time-self.__reftime)/(self.max_time-self.min_time)

        print('Fitting ephemeris position relative to star coordinate {}'.format(self.star.to_string('hmsdms')))
        self.ksi = np.polyfit(dt, da.to(u.km).value, 2)
        self.eta = np.polyfit(dt, dd.to(u.km).value, 2)
        print('ksi = aksi*t\u00b2 + bksi*t + cksi')
        print('eta = aeta*t\u00b2 + beta*t + ceta')
        print('t=(jd-{})/({}-{})'.format(self.__reftime.jd, self.max_time.jd,self.min_time.jd))
        print('        aksi={}'.format(self.ksi[0]))
        print('        bksi={}'.format(self.ksi[1]))
        print('        cksi={}'.format(self.ksi[2]))
        print('        aeta={}'.format(self.eta[0]))
        print('        beta={}'.format(self.eta[1]))
        print('        ceta={}'.format(self.eta[2]))
        ksi = np.poly1d(self.ksi)
        eta = np.poly1d(self.eta)
        dksi = da.to(u.km).value - ksi(dt)
        deta = dd.to(u.km).value - eta(dt)
        rmsk = np.sqrt(np.mean(np.square(dksi)))
        rmse = np.sqrt(np.mean(np.square(deta)))
        print('Residual RMS: ksi={:.3f} km, eta={:.3f} km'.format(rmsk, rmse))
              
        
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
        if not time.isscalar:
            if any(time < self.min_time) or any(time > self.max_time):
                raise ValueError('time must be in the interval [{},{}]'.format(self.min_time, self.max_time))
        elif time.isscalar:
            if time < self.min_time or time > self.max_time:
                raise ValueError('time must be in the interval [{},{}]'.format(self.min_time, self.max_time))
        if hasattr(self, 'ksi') and hasattr(self, 'eta'):
            ksi = np.poly1d(self.ksi)
            eta = np.poly1d(self.eta)
            return ksi((time-self.__reftime)/(self.max_time-self.min_time)), eta((time-self.__reftime)/(self.max_time-self.min_time))
        else:
            raise ValueError('A "star" parameter is missing. Please run fit_d2_ksi_eta first.')

    def __str__(self):
        """ String representation of the EphemPlanete Class.
        """
        out = 'Ephemeris of {}.\n'.format(self.name)
        out += 'Valid from {} until {}\n'.format(self.min_time.iso, self.max_time.iso)
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
        