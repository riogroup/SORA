from astropy.coordinates import SkyCoord, EarthLocation, GCRS
from astropy.coordinates.matrix_utilities import rotation_matrix
from astropy.time import Time
import astropy.units as u
from astroquery.mpc import MPC
import numpy as np
from .config import test_attr
from .star import Star

def search_code_mpc():
    obs = MPC.get_observatory_codes()
    observatories = {}
    for line in obs:
        code = line['Code']
        lon = line['Longitude']*u.deg
        rcphi = line['cos']*6378.137*u.km
        rsphi = line['sin']*6378.137*u.km
        name = line['Name']
        site = EarthLocation.from_geocentric(rcphi*np.cos(lon), rcphi*np.sin(lon), rsphi)
        observatories[code] = (name, site)
    return observatories
        

class Observer():
    '''
    Docstring
    Define the observer object
    '''
    __names = []
    def __init__(self, *args, **kwargs):
        # Run initial parameters
        if 'name' not in kwargs:
            kwargs['name'] = 'User'
        if len(args) == 1 or 'code' in kwargs:
            if len(args) == 1:
                code = args[0]
                try:
                    code = test_attr(args[0], float, 'code')
                    code = '{:03.0f}'.format(code)
                    kwargs['code'] = code
                except:
                    pass
            self.code = kwargs['code']
            try:
                self.__name, self.site = search_code_mpc()[self.code]
            except:
                raise ValueError('code {} could not be located in MPC database'.format(self.code))
        elif len(args) == 2 or all(i in kwargs for i in ['name', 'site']):
            if len(args) == 2:
                kwargs['name'] = args[0]
                kwargs['site'] = args[1]
            self.__name = kwargs['name']
            self.site = test_attr(kwargs['site'], EarthLocation, 'site')
        elif len(args) == 4 or all(i in kwargs for i in ['name', 'lon', 'lat', 'height']):
            if len(args) == 4:
                kwargs['name'] = args[0]
                kwargs['lon'] = args[1]
                kwargs['lat'] = args[2]
                kwargs['height'] = args[3]
            self.__name = kwargs['name']
            self.site = EarthLocation(kwargs['lon'], kwargs['lat'], kwargs['height'])
        else:
            raise ValueError('Input parameters could not be determined')
        if self.__name in self.__names:
            raise ValueError('name {} already defined for another Observer object. Please choose a different one.'.format(self.__name))
        self.__names.append(self.__name)
        
    def get_ksi_eta(self, time, star):
        """ Calculates relative position to star in the orthographic projection.
        
        Parameters:
        time (str, Time):Time from which to calculate the position.
        It can be a string in the format "yyyy-mm-dd hh:mm:ss.s" or an astropy Time object
        star (str, SkyCoord):The coordinate of the star in the same frame as the ephemeris.
        It can be a string in the format "hh mm ss.s +dd mm ss.ss" or an astropy SkyCoord object.
        
        Returns:
        ksi, eta (float): on-sky orthographic projection of the observer relative to a star
        """
        time = test_attr(time, Time, 'time')
        try:
            star = SkyCoord(star, unit=(u.hourangle,u.deg))
        except:
            raise ValueError('star is not an astropy object or a string in the format "hh mm ss.s +dd mm ss.ss"')

        itrs = self.site.get_itrs(obstime=time)
        gcrs = itrs.transform_to(GCRS(obstime=time))
        rz = rotation_matrix(star.ra, 'z')
        ry = rotation_matrix(-star.dec, 'y')

        cp = gcrs.cartesian.transform(rz).transform(ry)
        return cp.y.to(u.km).value, cp.z.to(u.km).value
    
    def sidereal_time(self, time, mode='local'):
        # return local or greenwich sidereal time
        time = test_attr(time, Time, 'time')
        time.location = self.site
        if mode == 'local':
            return time.sidereal_time('apparent')
        elif mode == 'greenwich':
            return time.sidereal_time('apparent', 'greenwich')
        else:
            raise ValueError('mode must be "local" or "greenwich"')

    @property
    def name(self):
        return self.__name

    @property
    def lon(self):
        return self.site.lon

    @lon.setter
    def lon(self, lon):
        lat = self.site.lat
        height = self.site.height
        site = EarthLocation(lon, lat, height)
        self.site = site

    @property
    def lat(self):
        return self.site.lat

    @lat.setter
    def lat(self, lat):
        lon = self.site.lon
        height = self.site.height
        site = EarthLocation(lon, lat, height)
        self.site = site

    @property
    def height(self):
        return self.site.height

    @height.setter
    def height(self, height):
        lon = self.site.lon
        lat = self.site.lat
        site = EarthLocation(lon, lat, height)
        self.site = site

    def __str__(self):
        """String representation of the Observer class
        """
        out = 'Site: {}\n'.format(self.name)
        out += 'Geodetic coordinates: Lon: {}, Lat: {}, height: {:.3f}'.format(
            self.site.lon.__str__(), self.site.lat.__str__(), self.site.height.to(u.km))
        return out

    def __del__(self):
        try:
            self.__names.remove(self.__name)
        except:
            pass
