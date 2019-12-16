from astropy.coordinates import EarthLocation, GCRS
from astropy.coordinates.matrix_utilities import rotation_matrix
from astropy.time import Time
import astropy.units as u
import numpy as np
from .config import test_attr
from .star import Star

def search_code_mpc():
    from urllib.request import urlopen
    data = urlopen('https://www.minorplanetcenter.net/iau/lists/ObsCodes.html')
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
    def __init__(self, **kwargs):
        # Run initial parameters
        #if 'site' in kwargs:
        #    self.site = test_attr(kwargs['site'], EarthLocation, 'site')
        #elif 'lat' in kwargs and 'lon' in kwargs:
        #    height = 0.0
        #    if 'height' in kwargs:
        #        height = kwargs['height']
        #    self.site = EarthLocation(kwargs['lon'], kwargs['lat'], height)
        #self.name = ''
        #if 'name' in kwargs:
        #    self.name = test_attr(kwargs['name'], str, 'name')
        return
        
    def get_parallax(self, time, star):
        # return relative position to star in the orthographic projection.
        #time = test_attr(time, Time, 'time')
        #star = test_attr(star, Star, 'star')
        #time.location = self.site

        #itrs = self.site.get_itrs(obstime=time)
        #gcrs = itrs.transform_to(GCRS(obstime=time))
        #rz = rotation_matrix((star.ra-self.sidereal_time(time, 'greenwich')), 'z')
        #ry = rotation_matrix(-star.dec, 'y')

        #cp = itrs.cartesian.transform(rz).transform(ry)
        #return cp.y.to(u.km).value, cp.z.to(u.km).value
        return
    
    def sidereal_time(self, time, mode='local'):
        # return local or greenwich sidereal time
        return
            
    def __str__(self):
        # return what it is to be printed
        return ''