from astropy.coordinates import SkyCoord, EarthLocation, GCRS
from astropy.coordinates.matrix_utilities import rotation_matrix
from astropy.time import Time
import astropy.units as u
import numpy as np
from .config import test_attr
from .star import Star

def search_code_mpc():
    from urllib.request import urlopen
    data = urlopen('https://www.minorplanetcenter.net/iau/lists/ObsCodes.html')
    print('Looking the MPC database ...')
    lines = data.readlines()
    del lines[0:2]
    del lines[-1]
    observatories = {}
    for line in lines:
        try:
            obs = line.decode('utf-8').strip()
            code = obs[0:3]
            lon = float(obs[3:13])*u.deg
            rcphi = float(obs[13:21])*u.R_earth
            rsphi = float(obs[21:30])*u.R_earth
            name = obs[30:]
            site = EarthLocation(rcphi*np.cos(lon), rcphi*np.sin(lon), rsphi)
            observatories[code] = (name, site)
        except:
            pass
    return observatories
        

class Observer():
    '''
    Docstring
    Define the observer object
    '''
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
                self.name, self.site = search_code_mpc()[self.code]
            except:
                raise ValueError('code {} could not be located in MPC database'.format(self.code))
        elif len(args) == 2 or all(i in kwargs for i in ['name', 'site']):
            if len(args) == 2:
                kwargs['name'] = args[0]
                kwargs['site'] = args[1]
            self.name = kwargs['name']
            self.site = test_attr(kwargs['site'], EarthLocation, 'site')
        elif len(args) == 4 or all(i in kwargs for i in ['name', 'lon', 'lat', 'height']):
            if len(args) == 4:
                kwargs['name'] = args[0]
                kwargs['lon'] = args[1]
                kwargs['lat'] = args[2]
                kwargs['height'] = args[3]
            self.name = kwargs['name']
            self.site = EarthLocation(kwargs['lon'], kwargs['lat'], kwargs['height'])
        else:
            raise ValueError('Input parameters could not be determined')
        
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

    def __str__(self):
        # return what it is to be printed
        out = ''
        return out