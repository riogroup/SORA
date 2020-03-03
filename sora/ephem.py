import numpy as np
from astropy.coordinates import SkyCoord, SkyOffsetFrame
from astropy.time import Time
import astropy.units as u
import astropy.constants as const
from astroquery.jplhorizons import Horizons
from .config import test_attr
import os
import spiceypy as spice

class EphemPlanete():
    """ EphemPlanete simulates ephem_planete and fit_d2_ksi_eta.

    Parameters:
        name (str): name of the object for search in the JPL database
        ephem (str):Input file with hour, minute, RA, DEC, distance  

    """
    def __init__(self, name, ephem, **kwargs):
        data = np.loadtxt(ephem, unpack=True)
        self.name = name
        self.time = Time(data[0], format='jd')
        self.ephem = SkyCoord(data[1]*u.deg, data[2]*u.deg, data[3]*u.AU)
        self.__reftime = self.time[0]
        self.min_time = Time(data[0].min(), format='jd')
        self.max_time = Time(data[0].max(), format='jd')
        self.radius = 0*u.km
        if 'radius' in kwargs:
            self.radius = kwargs['radius']*u.km
        self.mass = 0*u.kg
        if 'mass' in kwargs:
            self.mass = kwargs['mass']*u.kg

    def fit_d2_ksi_eta(self, star, log=True):
        """ Fits the on-sky ephemeris position relative to a star
        
        Parameters:
        star (str, SkyCoord):The coordinate of the star in the same frame as the ephemeris.
        log (bool): if fit log is to printed. Default: True
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

        
        self.ksi = np.polyfit(dt, da.to(u.km).value, 2)
        self.eta = np.polyfit(dt, dd.to(u.km).value, 2)
        ksi = np.poly1d(self.ksi)
        eta = np.poly1d(self.eta)
        dksi = da.to(u.km).value - ksi(dt)
        deta = dd.to(u.km).value - eta(dt)
        rmsk = np.sqrt(np.mean(np.square(dksi)))
        rmse = np.sqrt(np.mean(np.square(deta)))
        if log:
            print('Fitting ephemeris position relative to star coordinate {}'.format(self.star.to_string('hmsdms')))
            print('ksi = aksi*t\u00b2 + bksi*t + cksi')
            print('eta = aeta*t\u00b2 + beta*t + ceta')
            print('t=(jd-{})/({}-{})'.format(self.__reftime.jd, self.max_time.jd,self.min_time.jd))
            print('        aksi={}'.format(self.ksi[0]))
            print('        bksi={}'.format(self.ksi[1]))
            print('        cksi={}'.format(self.ksi[2]))
            print('        aeta={}'.format(self.eta[0]))
            print('        beta={}'.format(self.eta[1]))
            print('        ceta={}'.format(self.eta[2]))
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
    
    
class EphemJPL():
    """ EphemJPL obtains the ephemeris from Horizons website.

    Parameters:
        name (str): name of the object for search in the JPL database
        id_type (str): type of object options: 'smallbody', 'majorbody'
        (planets but also anything that is not a small body), 'designation',
        'name', 'asteroid_name', 'comet_name', 'id' (Horizons id number),
        or 'smallbody' (find the closest match under any id_type), default: 'smallbody'

    """
    def __init__(self, name, id_type='smallbody', **kwargs):
        self.name = name
        self.id_type='majorbody'
        self.radius = 0*u.km
        if 'radius' in kwargs:
            self.radius = kwargs['radius']*u.km
        self.mass = 0*u.kg
        if 'mass' in kwargs:
            self.mass = kwargs['mass']*u.kg

    def get_position(self, time):
        """ Returns the geocentric position of the object.

        Parameters:
        time (int, float):Time from which to calculate the position.

        Returns:
        coord (SkyCoord): Astropy SkyCoord object with the coordinate at given time
        """
        obj = Horizons(id=self.name, id_type=self.id_type, location='geo', epochs=time.jd)
        eph = obj.ephemerides(extra_precision=True)
        coord = SkyCoord(eph['RA'], eph['DEC'], eph['delta'], frame='icrs', obstime=time)
        if len(coord) == 1:
            return coord[0]
        return coord

    def get_ksi_eta(self, time, star=None):
        """ Returns the on-sky position of the ephemeris relative to a star.

        Parameters:
        time (int, float):Time from which to calculate the position.
        star (str, SkyCoord):The coordinate of the star in the same frame as the ephemeris.

        Returns:
        ksi, eta (float): on-sky position of the ephemeris relative to a star
        """
        if type(star) == str:
            star = SkyCoord(star, unit=(u.hourangle, u.deg))
        coord = self.get_position(time)
        target = coord.transform_to(SkyOffsetFrame(origin=star))
        da = -target.cartesian.y
        dd = -target.cartesian.z
        return da.to(u.km).value, dd.to(u.km).value

    def __str__(self):
        """ String representation of the EphemPlanete Class.
        """
        out = 'Ephemeris of {}.'.format(self.name)
        return out


class EphemKernel():
    """ EphemHorizons gets online the ephemeris for an object.

    Parameters:
        name (str): name of the object for search in the JPL database
        code (str): kernel code of the targeting object
        list of paths for kernels

    """
    def __init__(self, name, code, *args, **kwargs):
        self.name = name
        self.code = str(code)
        for arg in args:
            spice.furnsh(arg)
        #spice.spkpos(self.code, 0, 'J2000', 'NONE', '399')
        spice.kclear()
        self.__kernels = args
        self.radius = 0*u.km
        if 'radius' in kwargs:
            self.radius = kwargs['radius']*u.km
        self.mass = 0*u.kg
        if 'mass' in kwargs:
            self.mass = kwargs['mass']*u.kg

    
    def get_position(self, time):
        """ Returns the geocentric position of the object.

        Parameters:
        time (int, float):Time from which to calculate the position.

        Returns:
        coord (SkyCoord): Astropy SkyCoord object with the coordinate at given time
        """
        for arg in self.__kernels:
            spice.furnsh(arg)
        t0 = Time('J2000', scale='tdb')
        if time.isscalar:
            time = Time([time])
        dt = (time - t0)
        delt = 0*u.s
        position1 = np.array(spice.spkpos('0', dt.sec, 'J2000', 'NONE', '399')[0])
        while True:
            ### calculates new time
            tempo = dt - delt
            ### calculates vector Solar System Baricenter -> Object
            position2 = spice.spkpos(self.code, tempo.sec, 'J2000', 'NONE', '0')[0]
            position = (position1 + position2).T
            ### calculates linear distance Earth Topocenter -> Object
            dist = np.linalg.norm(position, axis=0)*u.km
            ### calculates new light time
            delt = (dist/const.c).decompose()
            ### if difference between new and previous light time is smaller than 0.001 sec, than continue.
            if all(np.absolute(((dt - tempo) - delt).sec) < 0.001):
                break
        coord = SkyCoord(position[0], position[1], position[2], frame='icrs', unit=u.km, representation_type='cartesian', obstime=time)
        spice.kclear()
        coord_rd = SkyCoord(ra=coord.spherical.lon, dec=coord.spherical.lat, distance=coord.spherical.distance, obstime=time)
        if len(coord) == 1:
            return coord_rd[0]
        return coord_rd
    
    def get_ksi_eta(self, time, star=None):
        """ Returns the on-sky position of the ephemeris relative to a star.

        Parameters:
        time (int, float):Time from which to calculate the position.
        star (str, SkyCoord):The coordinate of the star in the same frame as the ephemeris.

        Returns:
        ksi, eta (float): on-sky position of the ephemeris relative to a star
        """
        if type(star) == str:
            star = SkyCoord(star, unit=(u.hourangle, u.deg))
        coord = self.get_position(time)
        target = coord.transform_to(SkyOffsetFrame(origin=star))  
        da = -target.cartesian.y
        dd = -target.cartesian.z
        return da.to(u.km).value, dd.to(u.km).value
    
    def __str__(self):
        """ String representation of the EphemPlanete Class.
        """
        out = 'Ephemeris of {}.\n'.format(self.name)
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
        