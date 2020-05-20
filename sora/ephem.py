import numpy as np
from astropy.coordinates import SkyCoord, SkyOffsetFrame, SphericalCosLatDifferential, ICRS
from astropy.time import Time
from astropy.table import vstack
import astropy.units as u
import astropy.constants as const
from astroquery.jplhorizons import Horizons
from .config import test_attr
import os
import spiceypy as spice
import urllib.request
import warnings

def read_obj_data():
    """ Reads online table with physical parameter for selected objects
    at http://devel2.linea.gov.br/~altair.gomes/radius.txt

    Return:
        python dictionary
    """
    obj = {}
    try:
        data = urllib.request.urlopen('http://devel2.linea.gov.br/~altair.gomes/radius.txt')
        for line in data:
            arr = line.split()
            obj[arr[0].decode().lower()] = [float(i) for i in arr[1:]]
    except:
        warnings.warn('Online object data table could not be found. Please check internet connection.')
    return obj


def apparent_mag(H, G, dist, sundist, phase=0.0):
    """ Calculates the Apparent Magnitude

    Parameters:
        H (int, float): Absolute Magnitude
        G (int, float): Slope parameter
        dist (int, float): Observer-Object distance, in AU
        sundist (int, float): Sun-Object distance, in AU
        phase (int, float): Phase Angle: Sun-Target-Observer, in deg

    Return:
        ap_mag (float): Apparent Magnitude
    """
    phi0 = np.exp(-3.33*(np.tan(0.5*phase*u.deg)**0.63))
    phi1 = np.exp(-1.87*(np.tan(0.5*phase*u.deg)**1.22))
    Ha = H - 2.5*np.log10((1.0-G)*phi0 + G*phi1)
    ap_mag = Ha + 5*np.log10(sundist*dist)
    return ap_mag.value


def ephem_kernel(time, target, observer, kernels):
    """Calculate the ephemeris from kernel files

    Parameters:
        time (str, Time): instant to calculate ephemeris
        target (str): IAU (kernel) code of the target
        observer (str): IAU (kernel) code of the observer
        kernels (list, str): list of paths for all the kernels

    Return:
        coord (SkyCoord): ICRS coordinate of the target.
    """
    for kern in kernels:
        spice.furnsh(kern)
    time = Time(time)
    t0 = Time('J2000', scale='tdb')
    if time.isscalar:
        time = Time([time])
    dt = (time - t0)
    delt = 0*u.s
    position1 = np.array(spice.spkpos('0', dt.sec, 'J2000', 'NONE', observer)[0])
    while True:
        ### calculates new time
        tempo = dt - delt
        ### calculates vector Solar System Baricenter -> Object
        position2 = spice.spkpos(target, tempo.sec, 'J2000', 'NONE', '0')[0]
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


class EphemPlanete():
    """ EphemPlanete simulates ephem_planete and fit_d2_ksi_eta.

    Parameters:
        name (str): name of the object for search in the JPL database
        ephem (str):Input file with JD, RA, DEC, distance
        radius (int,float): Object radius, in km (Default: Online database)
        error_ra (int,float): Ephemeris RA*cosDEC error, in mas (Default: Online database)
        error_dec (int,float): Ephemeris DEC error, in mas (Default: Online database)
        mass (int,float): Object Mass, in kg (Default: 0.0)
        H (int,float): Object Absolute Magnitude (Default: NaN)
        G (int,float): Object Phase slope (Default: NaN)
    """
    def __init__(self, name, ephem, **kwargs):
        data = np.loadtxt(ephem, unpack=True)
        self.name = name
        self.time = Time(data[0], format='jd')
        self.ephem = SkyCoord(data[1]*u.deg, data[2]*u.deg, data[3]*u.AU)
        self.__reftime = self.time[0]
        self.min_time = Time(data[0].min(), format='jd')
        self.max_time = Time(data[0].max(), format='jd')
        try:
            data = read_obj_data()
        except:
            data = {}
        radius, error_ra, error_dec = data.get(name.lower(), [0,0,0])
        self.radius = kwargs.get('radius', radius)*u.km
        self.error_ra = kwargs.get('error_ra', error_ra)*u.arcsec
        self.error_dec = kwargs.get('error_dec', error_dec)*u.arcsec
        self.mass = kwargs.get('mass', 0.0)*u.kg
        self.H = kwargs.get('H', np.nan)
        self.G = kwargs.get('G', np.nan)

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
        da = target.cartesian.y
        dd = target.cartesian.z
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
        time = Time(time)
        if not time.isscalar:
            if any(time < self.min_time) or any(time > self.max_time):
                raise ValueError('time must be in the interval [{},{}]'.format(self.min_time, self.max_time))
        elif time.isscalar:
            if time < self.min_time or time > self.max_time:
                raise ValueError('time must be in the interval [{},{}]'.format(self.min_time, self.max_time))
        if hasattr(self, 'ksi') and hasattr(self, 'eta'):
            ksi = np.poly1d(self.ksi)
            eta = np.poly1d(self.eta)
            k = ksi((time-self.__reftime)/(self.max_time-self.min_time))
            e = eta((time-self.__reftime)/(self.max_time-self.min_time))
            if hasattr(self, 'offset'):
                dist = self.ephem[0].distance.to(u.km).value
                dksi = np.sin(self.offset.d_lon_coslat)*dist
                deta = np.sin(self.offset.d_lat)*dist
                k = k + dksi
                e = e + deta
            return k.value, e.value
        else:
            raise ValueError('A "star" parameter is missing. Please run fit_d2_ksi_eta first.')

    def apparent_magnitude(self,time):
        """ Calculates the Apparent Diameter

        Parameters:
            time (int, float):Time from which to calculate the magnitude.

        Returns:
            ap_mag (float): Apparent magnitude
        """
        time = Time(time)
        obj = Horizons(id=self.name, id_type='majorbody', location='geo', epochs=time.jd)
        eph = obj.ephemerides(extra_precision=True)
        if 'H' in eph.keys():
            self.H = eph['H'][0]
            self.G = eph['G'][0]
        if len(eph['V']) == 1:
            return eph['V'][0]
        else:
            return eph['V'].tolist()

    def add_offset(self, da_cosdec, ddec):
        """Add an offset to the Ephemeris

        Parameters:
            da_cosdec (int, float):Delta_alpha_cos_delta in mas
            ddec (int, float):Delta_delta in mas
        """
        dadc = test_attr(da_cosdec, float, 'da_cosdec')
        dd = test_attr(ddec, float, 'ddec')
        self.offset = SphericalCosLatDifferential(dadc*u.mas, dd*u.mas, 0.0*u.km)

    def __str__(self):
        """ String representation of the EphemPlanete Class.
        """
        out = 'Ephemeris of {}.\n'.format(self.name)
        out += 'Radius: {:.1f}\n'.format(self.radius)
        out += 'Mass: {:.2e}\n'.format(self.mass)
        out += '\nEphem Planete:\n'
        out += 'Valid from {} until {}\n'.format(self.min_time.iso, self.max_time.iso)
        out += 'Ephem Error: RA*cosDEC: {:.3f}; DEC: {:.3f}\n'.format(self.error_ra, self.error_dec)
        if hasattr(self,'star'):
            out += 'Fitted ephemeris position relative to star coordinate {}\n'.format(self.star.to_string('hmsdms'))
            out += '    ksi = aksi*t\u00b2 + bksi*t + cksi\n'
            out += '    eta = aeta*t\u00b2 + beta*t + ceta\n'
            out += '    t=(jd-{})/({}-{})\n'.format(self.__reftime.jd, self.max_time.jd,self.min_time.jd)
            out += '        aksi={}\n'.format(self.ksi[0])
            out += '        bksi={}\n'.format(self.ksi[1])
            out += '        cksi={}\n'.format(self.ksi[2])
            out += '        aeta={}\n'.format(self.eta[0])
            out += '        beta={}\n'.format(self.eta[1])
            out += '        ceta={}\n'.format(self.eta[2])
            target = self.ephem.transform_to(SkyOffsetFrame(origin=self.star))
            da = -target.cartesian.y
            dd = -target.cartesian.z
            dt = (self.time-self.__reftime)/(self.max_time-self.min_time)
            ksi = np.poly1d(self.ksi)
            eta = np.poly1d(self.eta)
            dksi = da.to(u.km).value - ksi(dt)
            deta = dd.to(u.km).value - eta(dt)
            rmsk = np.sqrt(np.mean(np.square(dksi)))
            rmse = np.sqrt(np.mean(np.square(deta)))
            out += 'Residual RMS: ksi={:.3f} km, eta={:.3f} km\n'.format(rmsk, rmse)
        if hasattr(self, 'offset'):
            out += '\nOffset applied: RA*cosDEC: {:.4f}; DEC: {:.4f}\n'.format(
                self.offset.d_lon_coslat.to(u.arcsec), self.offset.d_lat.to(u.arcsec))
        return out


class EphemJPL():
    """ EphemJPL obtains the ephemeris from Horizons service.

    Parameters:
        name (str): name of the object for search in the JPL database
        id_type (str): type of object options: 'smallbody', 'majorbody'
        (planets but also anything that is not a small body), 'designation',
        'name', 'asteroid_name', 'comet_name', 'id' (Horizons id number),
        or 'smallbody' (find the closest match under any id_type), default: 'smallbody'
        radius (int,float): Object radius, in km (Default: Online database)
        error_ra (int,float): Ephemeris RA*cosDEC error, in mas (Default: Online database)
        error_dec (int,float): Ephemeris DEC error, in mas (Default: Online database)
        mass (int,float): Object Mass, in kg (Default: 0.0)
        H (int,float): Object Absolute Magnitude (Default: NaN)
        G (int,float): Object Phase slope (Default: NaN)
    """
    def __init__(self, name, id_type='smallbody', **kwargs):
        self.name = name
        self.id_type='majorbody'
        try:
            data = read_obj_data()
        except:
            data = {}
        radius, error_ra, error_dec = data.get(name.lower(), [0,0,0])
        self.radius = kwargs.get('radius', radius)*u.km
        self.error_ra = kwargs.get('error_ra', error_ra)*u.arcsec
        self.error_dec = kwargs.get('error_dec', error_dec)*u.arcsec
        self.mass = kwargs.get('mass', 0.0)*u.kg
        self.H = kwargs.get('H', np.nan)
        self.G = kwargs.get('G', np.nan)

    def get_position(self, time):
        """ Returns the geocentric position of the object.

        Parameters:
        time (int, float):Time from which to calculate the position.

        Returns:
        coord (SkyCoord): Astropy SkyCoord object with the coordinate at given time
        """
        time = Time(time)
        if not time.isscalar:
            s = len(time)
            def calc_ephem(i):
                n = 50*i
                k = n+50
                if k > s:
                    k = s
                obj = Horizons(id=self.name, id_type=self.id_type, location='geo', epochs=time[n:k].jd)
                return obj.ephemerides(extra_precision=True)
            plus = 0
            if s%50 > 0:
                plus = 1
            eph = vstack([calc_ephem(j) for j in range(s//50+1*plus)])
        else:
            obj = Horizons(id=self.name, id_type=self.id_type, location='geo', epochs=time.jd)
            eph = obj.ephemerides(extra_precision=True)
        coord = SkyCoord(eph['RA'], eph['DEC'], eph['delta'], frame='icrs', obstime=time)
        if hasattr(self, 'offset'):
            pos_frame = SkyOffsetFrame(origin=coord)
            new_pos = SkyCoord(lon=self.offset.d_lon_coslat, lat=self.offset.d_lat, distance=coord.distance, frame=pos_frame)
            coord = new_pos.transform_to(ICRS)
        if len(coord) == 1:
            return coord[0]
        return coord

    def apparent_magnitude(self,time):
        """ Calculates the Apparent Diameter

        Parameters:
            time (int, float):Time from which to calculate the magnitude.

        Returns:
            ap_mag (float): Apparent magnitude
        """
        time = Time(time)
        obj = Horizons(id=self.name, id_type='majorbody', location='geo', epochs=time.jd)
        eph = obj.ephemerides(extra_precision=True)
        if 'H' in eph.keys():
            self.H = eph['H'][0]
            self.G = eph['G'][0]
        if len(eph['V']) == 1:
            return eph['V'][0]
        else:
            return eph['V'].tolist()

    def get_ksi_eta(self, time, star):
        """ Returns the on-sky position of the ephemeris relative to a star.

        Parameters:
        time (int, float):Time from which to calculate the position.
        star (str, SkyCoord):The coordinate of the star in the same frame as the ephemeris.

        Returns:
        ksi, eta (float): on-sky position of the ephemeris relative to a star
        """
        time = Time(time)
        if type(star) == str:
            star = SkyCoord(star, unit=(u.hourangle, u.deg))
        coord = self.get_position(time)
        target = coord.transform_to(SkyOffsetFrame(origin=star))
        da = target.cartesian.y
        dd = target.cartesian.z
        return da.to(u.km).value, dd.to(u.km).value

    def  get_pole_position_angle(self,pole,time):
        """ Returns the geocentric position of the object.

        Parameters:
        pole (str, astropy.SkyCoord):Coordinate of the pole ICRS.
        time (int, float): Time from which to calculate the position.

        Returns:
        position_angle (float): Position angle of the pole, in degrees
        aperture_angle (float): Apeture angle of the pole, in degrees
        """
        time = Time(time)
        if type(pole) == str:
            pole = SkyCoord(pole, unit=(u.hourangle, u.deg))
        obj = self.get_position(time)
        position_angle = obj.position_angle(pole).value*u.rad
        aperture_angle = np.arcsin(-1*(np.sin(pole.dec)*np.sin(obj.dec) + np.cos(pole.dec)*np.cos(obj.dec)*np.cos(pole.ra-obj.ra)))
        return position_angle.to('deg'), aperture_angle.to('deg')

    def add_offset(self, da_cosdec, ddec):
        """Add an offset to the Ephemeris

        Parameters:
            da_cosdec (int, float):Delta_alpha_cos_delta in mas
            ddec (int, float):Delta_delta in mas
        """
        dadc = test_attr(da_cosdec, float, 'da_cosdec')
        dd = test_attr(ddec, float, 'ddec')
        self.offset = SphericalCosLatDifferential(dadc*u.mas, dd*u.mas, 0.0*u.km)

    def __str__(self):
        """ String representation of the EphemPlanete Class.
        """
        out = 'Ephemeris of {}.\n'.format(self.name)
        out += 'Radius: {:.1f}\n'.format(self.radius)
        out += 'Mass: {:.2e}\n'.format(self.mass)
        out += '\nEphemeris are downloaded from Horizons website\n'
        if hasattr(self, 'offset'):
            out += 'Offset applied: RA*cosDEC: {:.4f}; DEC: {:.4f}\n'.format(
                self.offset.d_lon_coslat.to(u.arcsec), self.offset.d_lat.to(u.arcsec))
        return out


class EphemKernel():
    """ EphemKernel gets the ephemeris from bsp kernels.

    Parameters:
        name (str): name of the object for search in the JPL database
        code (str): kernel code of the targeting object
        kernels(list): list of paths for kernels
        radius (int,float): Object radius, in km (Default: Online database)
        error_ra (int,float): Ephemeris RA*cosDEC error, in mas (Default: Online database)
        error_dec (int,float): Ephemeris DEC error, in mas (Default: Online database)
        mass (int,float): Object Mass, in kg (Default: 0.0)
        H (int,float): Object Absolute Magnitude (Default: NaN)
        G (int,float): Object Phase slope (Default: NaN)
    """
    def __init__(self, name, code, kernels, **kwargs):
        self.name = name
        self.code = str(code)
        self.meta = {}
        kerns = []
        for arg in kernels:
            spice.furnsh(arg)
            kerns.append(arg.split('/')[-1].split('.')[0].upper())
        spice.kclear()
        self.meta['kernels'] = '/'.join(kerns)
        self.__kernels = kernels
        try:
            data = read_obj_data()
        except:
            data = {}
        radius, error_ra, error_dec = data.get(name.lower(), [0,0,0])
        self.radius = kwargs.get('radius', radius)*u.km
        self.error_ra = kwargs.get('error_ra', error_ra)*u.arcsec
        self.error_dec = kwargs.get('error_dec', error_dec)*u.arcsec
        self.mass = kwargs.get('mass', 0.0)*u.kg
        self.H = kwargs.get('H', np.nan)
        self.G = kwargs.get('G', np.nan)

    def get_position(self, time):
        """ Returns the geocentric position of the object.

        Parameters:
        time (int, float):Time from which to calculate the position.

        Returns:
        coord (SkyCoord): Astropy SkyCoord object with the coordinate at given time
        """
        pos = ephem_kernel(time, self.code, '399', self.__kernels)
        if hasattr(self, 'offset'):
            pos_frame = SkyOffsetFrame(origin=pos)
            new_pos = SkyCoord(lon=self.offset.d_lon_coslat, lat=self.offset.d_lat, distance=pos.distance, frame=pos_frame)
            pos = new_pos.transform_to(ICRS)
        return pos

    def apparent_magnitude(self,time):
        """ Calculates the Apparent Diameter

        Parameters:
            time (int, float):Time from which to calculate the magnitude.

        Returns:
            ap_mag (float): Apparent magnitude
        """
        time = Time(time)

        if np.isnan(self.H) or np.isnan(self.G):
            obj = Horizons(id=self.name, id_type='majorbody', location='geo', epochs=time.jd)
            eph = obj.ephemerides(extra_precision=True)
            if 'H' in eph.keys():
                self.H = eph['H'][0]
                self.G = eph['G'][0]
            if len(eph['V']) == 1:
                return eph['V'][0]
            else:
                return eph['V'].tolist()

        else:
            obs_obj = self.get_position(time)
            sun_obj = ephem_kernel(time, self.code, '10', self.__kernels)

            # Calculates the phase angle between the 2-vectors
            unit_vector_1 = -obs_obj.cartesian.xyz / np.linalg.norm(obs_obj.cartesian.xyz)
            unit_vector_2 = -sun_obj.cartesian.xyz / np.linalg.norm(sun_obj.cartesian.xyz)
            dot_product = np.dot(unit_vector_1, unit_vector_2)
            phase = np.arccos(dot_product).to(u.deg).value

            return apparent_mag(self.H, self.G, obs_obj.distance.to(u.AU).value, sun_obj.distance.to(u.AU).value, phase)

    def get_ksi_eta(self, time, star):
        """ Returns the on-sky position of the ephemeris relative to a star.

        Parameters:
        time (int, float):Time from which to calculate the position.
        star (str, SkyCoord):The coordinate of the star in the same frame as the ephemeris.

        Returns:
        ksi, eta (float): on-sky position of the ephemeris relative to a star
        """
        time = Time(time)
        if type(star) == str:
            star = SkyCoord(star, unit=(u.hourangle, u.deg))
        coord = self.get_position(time)
        target = coord.transform_to(SkyOffsetFrame(origin=star))  
        da = target.cartesian.y
        dd = target.cartesian.z
        return da.to(u.km).value, dd.to(u.km).value

    def  get_pole_position_angle(self,pole,time):
        """ Returns the geocentric position of the object.

        Parameters:
        pole (str, astropy.SkyCoord):Coordinate of the pole ICRS.
        time (int, float): Time from which to calculate the position.

        Returns:
        position_angle (float): Position angle of the pole, in degrees
        aperture_angle (float): Apeture angle of the pole, in degrees
        """
        time = Time(time)
        if type(pole) == str:
            pole = SkyCoord(pole, unit=(u.hourangle, u.deg))
        obj = self.get_position(time)
        position_angle = obj.position_angle(pole).value*u.rad
        aperture_angle = np.arcsin(-1*(np.sin(pole.dec)*np.sin(obj.dec) + np.cos(pole.dec)*np.cos(obj.dec)*np.cos(pole.ra-obj.ra)))
        return position_angle.to('deg'), aperture_angle.to('deg')

    def add_offset(self, da_cosdec, ddec):
        """Add an offset to the Ephemeris

        Parameters:
            da_cosdec (int, float):Delta_alpha_cos_delta in mas
            ddec (int, float):Delta_delta in mas
        """
        dadc = test_attr(da_cosdec, float, 'da_cosdec')
        dd = test_attr(ddec, float, 'ddec')
        self.offset = SphericalCosLatDifferential(dadc*u.mas, dd*u.mas, 0.0*u.km)

    def __str__(self):
        """ String representation of the EphemPlanete Class.
        """
        out = 'Ephemeris of {}.\n'.format(self.name)
        out += 'Radius: {:.1f}\n'.format(self.radius)
        out += 'Mass: {:.2e}\n'.format(self.mass)
        out += '\nEphem Kernel: {} (code={})\n'.format(self.meta['kernels'],self.code)
        out += 'Ephem Error: RA*cosDEC: {:.3f}; DEC: {:.3f}\n'.format(self.error_ra, self.error_dec)
        if hasattr(self, 'offset'):
            out += 'Offset applied: RA*cosDEC: {:.4f}; DEC: {:.4f}\n'.format(
                self.offset.d_lon_coslat.to(u.arcsec), self.offset.d_lat.to(u.arcsec))
        return out
