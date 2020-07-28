import numpy as np
from astropy.coordinates import SkyCoord, SkyOffsetFrame, \
    SphericalCosLatDifferential, ICRS
from astropy.time import Time
from astropy.table import vstack
import astropy.units as u
import astropy.constants as const
from astroquery.jplhorizons import Horizons
from sora.config import test_attr, input_tests
from sora.config.decorators import deprecated_alias
import spiceypy as spice
import urllib.request
import warnings


warnings.simplefilter('always', UserWarning)


def read_obj_data():
    """ Reads an online table (link below) with physical parameters for selected objects

    Table url: http://devel2.linea.gov.br/~altair.gomes/radius.txt

    Table content: Object name; radius (km); uncertainty in RA; uncertainty in DEC
        RA: Delta * alpha * cos (delta)
        DEC: Delta * delta

    Returns:
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
    """ Calculates the Object Apparent Magnitude

    Parameters:
        H (int, float): Absolute Magnitude
        G (int, float): Slope parameter
        dist (int, float): Observer-Object distance, in AU
        sundist (int, float): Sun-Object distance, in AU
        phase (int, float): Phase Angle: Sun-Target-Observer, in deg

    Returns:
        ap_mag (float): Apparent Magnitude
    """
    phi0 = np.exp(-3.33*(np.tan(0.5*phase*u.deg)**0.63))
    phi1 = np.exp(-1.87*(np.tan(0.5*phase*u.deg)**1.22))
    Ha = H - 2.5*np.log10((1.0-G)*phi0 + G*phi1)
    ap_mag = Ha + 5*np.log10(sundist*dist)
    return ap_mag.value


def ephem_kernel(time, target, observer, kernels):
    """ Calculates the ephemeris from kernel files

    Parameters:
        time (str, Time): reference instant to calculate ephemeris
        target (str): IAU (kernel) code of the target
        observer (str): IAU (kernel) code of the observer
        kernels (list, str): list of paths for all the kernels

    Returns:
        coord (SkyCoord): ICRS coordinate of the target.
    """
    if type(kernels) == str:
        kernels = [kernels]
    for kern in kernels:
        spice.furnsh(kern)
    time = Time(time)
    t0 = Time('J2000', scale='tdb')
    if time.isscalar:
        time = Time([time])
    dt = (time - t0)
    delt = 0*u.s
    # calculates vector Observer -> Solar System Baricenter
    position1 = np.array(spice.spkpos('0', dt.sec, 'J2000', 'NONE', observer)[0])
    while True:
        # calculates new time
        tempo = dt - delt
        # calculates vector Solar System Baricenter -> Object
        position2 = spice.spkpos(target, tempo.sec, 'J2000', 'NONE', '0')[0]
        position = (position1 + position2).T
        # calculates linear distance Earth Topocenter -> Object
        dist = np.linalg.norm(position, axis=0)*u.km
        # calculates new light time
        delt = (dist/const.c).decompose()
        # if difference between new and previous light time is smaller than 0.001 sec, than continue.
        if all(np.absolute(((dt - tempo) - delt).sec) < 0.001):
            break
    coord = SkyCoord(position[0], position[1], position[2], frame='icrs', unit=u.km,
                     representation_type='cartesian', obstime=time)
    spice.kclear()
    coord_rd = SkyCoord(ra=coord.spherical.lon, dec=coord.spherical.lat,
                        distance=coord.spherical.distance, obstime=time)
    if len(coord) == 1:
        return coord_rd[0]
    return coord_rd


class EphemPlanete():
    def __init__(self, name, ephem, **kwargs):
        """ Simulates the former fortran programs ephem_planete and fit_d2_ksi_eta.

        Parameters:
            name (str): name of the object to search in the JPL database
            ephem (str): Input file with JD (UTC), and geocentric RA (deg), DEC (deg), and distance (AU)
            radius (int,float): Object radius, in km (Default: Online database)
            error_ra (int,float): Ephemeris RA*cosDEC error, in arcsec (Default: Online database)
            error_dec (int,float): Ephemeris DEC error, in arcsec (Default: Online database)
            mass (int,float): Object Mass, in kg (Default: 0.0)
            H (int,float): Object Absolute Magnitude (Default: NaN)
            G (int,float): Object Phase slope (Default: NaN)
        """
        input_tests.check_kwargs(kwargs, allowed_kwargs=['error_dec', 'error_ra', 'H', 'G', 'mass', 'radius'])
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
        radius, error_ra, error_dec = data.get(name.lower(), [0, 0, 0])
        self.radius = kwargs.get('radius', radius)*u.km
        self.error_ra = kwargs.get('error_ra', error_ra)*u.arcsec
        self.error_dec = kwargs.get('error_dec', error_dec)*u.arcsec
        self.mass = kwargs.get('mass', 0.0)*u.kg
        self.H = kwargs.get('H', np.nan)
        self.G = kwargs.get('G', np.nan)

    def get_position(self, time):
        """ Returns the geocentric position of the object.

        Parameters:
            time (str, Time):Time from which to calculate the position.

        Returns:
            coord (SkyCoord): Astropy SkyCoord object with the coordinate at given time
        """
        ksi, eta = self.get_ksi_eta(time=time)*u.km
        distance = self.ephem.distance.mean()
        off_ra = np.arctan2(ksi, distance)
        off_dec = np.arctan2(eta, distance)
        coord_frame = SkyOffsetFrame(origin=self.star)
        pos = SkyCoord(lon=off_ra, lat=off_dec, distance=distance, frame=coord_frame)
        return pos.icrs

    def fit_d2_ksi_eta(self, star, log=True):
        """ Fits the projected position* of the object in the tangent sky plane relative to a star
            * ortographic projection.

        Parameters:
            star (str, SkyCoord): The coordinate of the star in the same reference frame as the ephemeris.
            log (bool): if True, log is printed. Default: True
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
            output = ('Fitting ephemeris position relative to star coordinate {}\n'
                      'ksi = aksi*t\u00b2 + bksi*t + cksi\n'
                      'eta = aeta*t\u00b2 + beta*t + ceta\n'
                      't=(jd-{})/({}-{})\n'
                      '        aksi={}\n'
                      '        bksi={}\n'
                      '        cksi={}\n'
                      '        aeta={}\n'
                      '        beta={}\n'
                      '        ceta={}\n'
                      'Residual RMS: ksi={:.3f} km, eta={:.3f} km'.format(
                          self.star.to_string('hmsdms'), self.__reftime.jd, self.max_time.jd,
                          self.min_time.jd, *self.ksi, *self.eta, rmsk, rmse)
                      )
            print(output)

    def get_ksi_eta(self, time, star=None):
        """ Returns the projected position* of the object in the tangent sky plane relative to a star.
            * ortographic projection.

        Parameters:
            time (str, Time): Reference time to calculate the position.
            star (str, SkyCoord): The coordinate of the star in the same reference frame as the ephemeris.

        Returns:
            ksi, eta (float): projected position (ortographic projection) of the object in the tangent sky plane
                relative to a star.
                Ksi is in the North-South direction (North positive)
                Eta is in the East-West direction (East positive)
        """
        if star:
            self.fit_d2_ksi_eta(star)
        time = Time(time)
        if not time.isscalar:
            if any(time < self.min_time) or any(time > self.max_time):
                raise ValueError('time must be in the interval [{},{}]'.format(
                    self.min_time, self.max_time))
        elif time.isscalar:
            if time < self.min_time or time > self.max_time:
                raise ValueError('time must be in the interval [{},{}]'.format(
                    self.min_time, self.max_time))
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

    def apparent_magnitude(self, time):
        """ Calculates the object Apparent Magnitude

        Parameters:
            time (str, Time): reference time to calculate the apparent magnitude.

        Returns:
            ap_mag (float): Object apparent magnitude
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
        """ Adds an offset to the Ephemeris

        Parameters:
            da_cosdec (int, float): Delta_alpha_cos_delta in mas
            ddec (int, float): Delta_delta in mas
        """
        dadc = test_attr(da_cosdec, float, 'da_cosdec')
        dd = test_attr(ddec, float, 'ddec')
        self.offset = SphericalCosLatDifferential(dadc*u.mas, dd*u.mas, 0.0*u.km)

    def __str__(self):
        """ String representation of the EphemPlanete Class.
        """
        out = ('Ephemeris of {}.\n'
               'Radius: {:.1f}\n'
               'Mass: {:.2e}\n'
               '\nEphem Planete:\n'
               'Valid from {} until {}\n'
               'Ephem Error: RA*cosDEC: {:.3f}; DEC: {:.3f}\n'.format(
                   self.name, self.radius, self.mass, self.min_time.iso,
                   self.max_time.iso, self.error_ra, self.error_dec)
               )
        if hasattr(self, 'star'):
            out += ("Fitted ephemeris position relative to star coordinate {}\n"
                    "    ksi = aksi*t\u00b2 + bksi*t + cksi\n"
                    "    eta = aeta*t\u00b2 + beta*t + ceta\n"
                    "    t=(jd-{})/({}-{})\n"
                    "        aksi={}\n"
                    "        bksi={}\n"
                    "        aeta={}\n"
                    "        beta={}\n"
                    "        ceta={}\n".format(
                        self.star.to_string('hmsdms'), self.__reftime.jd, self.max_time.jd,
                        self.min_time.jd, *self.ksi, *self.eta)
                    )
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
    def __init__(self, name, id_type='smallbody', **kwargs):
        """ Obtains the ephemeris from Horizons/JPL service.

        Web tool URL = https://ssd.jpl.nasa.gov/horizons.cgi

        Parameters:
            name (str): name of the object to search in the JPL database
            id_type (str): type of object options: 'smallbody', 'majorbody'
                (planets but also anything that is not a small body), 'designation',
                'name', 'asteroid_name', 'comet_name', 'id' (Horizons id number),
                or 'smallbody' (find the closest match under any id_type). Default: 'smallbody'
            radius (int,float): Object radius, in km (Default: Online database)
            error_ra (int,float): Ephemeris RA*cosDEC error, in arcsec (Default: Online database)
            error_dec (int,float): Ephemeris DEC error, in arcsec (Default: Online database)
            mass (int,float): Object Mass, in kg (Default: 0.0)
            H (int,float): Object Absolute Magnitude (Default: NaN)
            G (int,float): Object Phase slope (Default: NaN)
        """
        input_tests.check_kwargs(kwargs, allowed_kwargs=['error_dec', 'error_ra', 'H', 'G', 'mass', 'radius'])
        self.name = name
        self.id_type = id_type
        try:
            data = read_obj_data()
        except:
            data = {}
        radius, error_ra, error_dec = data.get(name.lower(), [0, 0, 0])
        self.radius = kwargs.get('radius', radius)*u.km
        self.error_ra = kwargs.get('error_ra', error_ra)*u.arcsec
        self.error_dec = kwargs.get('error_dec', error_dec)*u.arcsec
        self.mass = kwargs.get('mass', 0.0)*u.kg
        self.H = kwargs.get('H', np.nan)
        self.G = kwargs.get('G', np.nan)

    def get_position(self, time):
        """ Returns the geocentric position of the object.

        Parameters:
            time (str, Time): Reference time to calculate the position.

        Returns:
            coord (SkyCoord): Astropy SkyCoord object with the object coordinates at the given time
        """
        time = Time(time)
        if not time.isscalar:
            s = len(time)

            def calc_ephem(i):
                n = 50*i
                k = n+50
                if k > s:
                    k = s
                obj = Horizons(id=self.name, id_type=self.id_type,
                               location='geo', epochs=time[n:k].jd)
                return obj.ephemerides(extra_precision=True)

            plus = 0
            if s % 50 > 0:
                plus = 1
            eph = vstack([calc_ephem(j) for j in range(s//50+1*plus)])
        else:
            obj = Horizons(id=self.name, id_type=self.id_type, location='geo', epochs=time.jd)
            eph = obj.ephemerides(extra_precision=True)
        coord = SkyCoord(eph['RA'], eph['DEC'], eph['delta'], frame='icrs', obstime=time)
        if hasattr(self, 'offset'):
            pos_frame = SkyOffsetFrame(origin=coord)
            new_pos = SkyCoord(lon=self.offset.d_lon_coslat, lat=self.offset.d_lat,
                               distance=coord.distance, frame=pos_frame)
            coord = new_pos.transform_to(ICRS)
        if len(coord) == 1:
            return coord[0]
        return coord

    def apparent_magnitude(self, time):
        """ Calculates the Object Apparent Magnitude

        Parameters:
            time (str, Time): Reference time to calculate the magnitude.

        Returns:
            ap_mag (float): Apparent magnitude
        """
        time = Time(time)
        obj = Horizons(id=self.name, id_type=self.id_type, location='geo', epochs=time.jd)
        eph = obj.ephemerides(extra_precision=True)
        if 'H' in eph.keys():
            self.H = eph['H'][0]
            self.G = eph['G'][0]
        if len(eph['V']) == 1:
            return eph['V'][0]
        else:
            return eph['V'].tolist()

    def get_ksi_eta(self, time, star):
        """ Returns projected position* of the object in the tangent sky plane relative to a star.
            * ortographic projection.


        Parameters:
            time (str, Time): Reference time to calculate the object position.
            star (str, SkyCoord): The coordinate of the star in the same reference frame as the ephemeris.

        Returns:
            ksi, eta (float): projected position (ortographic projection) of the object in the tangent sky plane
                relative to a star.
                Ksi is in the North-South direction (North positive)
                Eta is in the East-West direction (East positive)
        """
        time = Time(time)
        if type(star) == str:
            star = SkyCoord(star, unit=(u.hourangle, u.deg))
        coord = self.get_position(time)
        target = coord.transform_to(SkyOffsetFrame(origin=star))
        da = target.cartesian.y
        dd = target.cartesian.z
        return da.to(u.km).value, dd.to(u.km).value

    def get_pole_position_angle(self, pole, time):
        """ Returns the object geocentric position.

        Parameters:
            pole (str, SkyCoord): Coordinate of the object pole ICRS.
            time (str, Time): Reference time to calculate the position.

        Returns:
            position_angle (float): Position angle of the object pole, in degrees
            aperture_angle (float): Apeture angle of the object pole, in degrees
        """
        time = Time(time)
        if type(pole) == str:
            pole = SkyCoord(pole, unit=(u.hourangle, u.deg))
        obj = self.get_position(time)
        position_angle = obj.position_angle(pole).value*u.rad
        aperture_angle = np.arcsin(
            -(np.sin(pole.dec)*np.sin(obj.dec) +
              np.cos(pole.dec)*np.cos(obj.dec)*np.cos(pole.ra-obj.ra))
            )
        return position_angle.to('deg'), aperture_angle.to('deg')

    def add_offset(self, da_cosdec, ddec):
        """ Adds an offset to the Ephemeris

        Parameters:
            da_cosdec (int, float): Delta_alpha_cos_delta in mas
            ddec (int, float): Delta_delta in mas
        """
        dadc = test_attr(da_cosdec, float, 'da_cosdec')
        dd = test_attr(ddec, float, 'ddec')
        self.offset = SphericalCosLatDifferential(dadc*u.mas, dd*u.mas, 0.0*u.km)

    def __str__(self):
        """ String representation of the EphemJPL Class.
        """
        out = ("Ephemeris of {}.\n"
               "Radius: {:.1f}\n"
               "Mass: {:.2e}\n"
               "\nEphemeris are downloaded from Horizons website\n".format(
                   self.name, self.radius, self.mass)
               )
        if hasattr(self, 'offset'):
            out += 'Offset applied: RA*cosDEC: {:.4f}; DEC: {:.4f}\n'.format(
                self.offset.d_lon_coslat.to(u.arcsec), self.offset.d_lat.to(u.arcsec))
        return out


class EphemKernel():
    @deprecated_alias(code='spkid')  # remove this line for v1.0
    def __init__(self, name, spkid, kernels, **kwargs):
        """ Gets the ephemeris from bsp kernels.

        Parameters:
            name (str): name of the object to search in the JPL database
            spkid (str): spkid of the targeting object. Former 'code' (v0.1)
            kernels(list): list of paths for kernels files
            radius (int,float): Object radius, in km (Default: Online database)
            error_ra (int,float): Ephemeris RA*cosDEC error, in arcsec (Default: Online database)
            error_dec (int,float): Ephemeris DEC error, in arcsec (Default: Online database)
            mass (int,float): Object Mass, in kg (Default: 0.0)
            H (int,float): Object Absolute Magnitude (Default: NaN)
            G (int,float): Object Phase slope (Default: NaN)
        """
        input_tests.check_kwargs(kwargs, allowed_kwargs=['error_dec', 'error_ra', 'H', 'G', 'mass', 'radius'])
        self.name = name
        self.spkid = str(spkid)
        self.code = self.spkid  # remove this line for v1.0
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
        radius, error_ra, error_dec = data.get(name.lower(), [0, 0, 0])
        self.radius = kwargs.get('radius', radius)*u.km
        self.error_ra = kwargs.get('error_ra', error_ra)*u.arcsec
        self.error_dec = kwargs.get('error_dec', error_dec)*u.arcsec
        self.mass = kwargs.get('mass', 0.0)*u.kg
        self.H = kwargs.get('H', np.nan)
        self.G = kwargs.get('G', np.nan)

    def get_position(self, time):
        """ Returns the object geocentric position.

        Parameters:
            time (str, Time): Reference time to calculate the object position.

        Returns:
            coord (SkyCoord): Astropy SkyCoord object with the object coordinates at the given time
        """
        pos = ephem_kernel(time, self.spkid, '399', self.__kernels)
        if hasattr(self, 'offset'):
            pos_frame = SkyOffsetFrame(origin=pos)
            new_pos = SkyCoord(lon=self.offset.d_lon_coslat, lat=self.offset.d_lat,
                               distance=pos.distance, frame=pos_frame)
            pos = new_pos.transform_to(ICRS)
        return pos

    def apparent_magnitude(self, time):
        """ Calculates the Object Apparent Magnitude

        Parameters:
            time (str, Time): Reference time to calculate the object aparent magnitude.

        Returns:
            ap_mag (float): Object apparent magnitude
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
            sun_obj = ephem_kernel(time, self.spkid, '10', self.__kernels)

            # Calculates the phase angle between the 2-vectors
            unit_vector_1 = -obs_obj.cartesian.xyz / np.linalg.norm(obs_obj.cartesian.xyz)
            unit_vector_2 = -sun_obj.cartesian.xyz / np.linalg.norm(sun_obj.cartesian.xyz)
            dot_product = np.dot(unit_vector_1, unit_vector_2)
            phase = np.arccos(dot_product).to(u.deg).value

            return apparent_mag(self.H, self.G, obs_obj.distance.to(u.AU).value,
                                sun_obj.distance.to(u.AU).value, phase)

    def get_ksi_eta(self, time, star):
        """ Returns projected position* of the object in the tangent sky plane relative to a star.
            * ortographic projection.

        Parameters:
            time (str, Time): Reference time to calculate the object position.
            star (str, SkyCoord): Coordinate of the star in the same reference frame as the ephemeris.

        Returns:
            ksi, eta (float): projected position (ortographic projection) of the object in the tangent sky plane
                relative to a star.
                Ksi is in the North-South direction (North positive)
                Eta is in the East-West direction (East positive)
        """
        time = Time(time)
        if type(star) == str:
            star = SkyCoord(star, unit=(u.hourangle, u.deg))
        coord = self.get_position(time)
        target = coord.transform_to(SkyOffsetFrame(origin=star))
        da = target.cartesian.y
        dd = target.cartesian.z
        return da.to(u.km).value, dd.to(u.km).value

    def get_pole_position_angle(self, pole, time):
        """ Returns the object geocentric position.

        Parameters:
            pole (str, SkyCoord): Coordinate of the object pole ICRS.
            time (str, Time): Time from which to calculate the position.

        Returns:
            position_angle (float): Position angle of the object pole, in degrees
            aperture_angle (float): Apeture angle of the object pole, in degrees
        """
        time = Time(time)
        if type(pole) == str:
            pole = SkyCoord(pole, unit=(u.hourangle, u.deg))
        obj = self.get_position(time)
        position_angle = obj.position_angle(pole).value*u.rad
        aperture_angle = np.arcsin(
            -(np.sin(pole.dec)*np.sin(obj.dec) +
              np.cos(pole.dec)*np.cos(obj.dec)*np.cos(pole.ra-obj.ra))
            )
        return position_angle.to('deg'), aperture_angle.to('deg')

    def add_offset(self, da_cosdec, ddec):
        """ Adds an offset to the Ephemeris

        Parameters:
            da_cosdec (int, float): Delta_alpha_cos_delta in mas
            ddec (int, float): Delta_delta in mas
        """
        dadc = test_attr(da_cosdec, float, 'da_cosdec')
        dd = test_attr(ddec, float, 'ddec')
        self.offset = SphericalCosLatDifferential(dadc*u.mas, dd*u.mas, 0.0*u.km)

    def __str__(self):
        """ String representation of the EphemKernel Class.
        """
        out = ("Ephemeris of {}.\n"
               "Radius: {:.1f}\n"
               "Mass: {:.2e}\n"
               "\nEphem Kernel: {} (SPKID={})\n"
               "Ephem Error: RA*cosDEC: {:.3f}; DEC: {:.3f}\n".format(
                   self.name, self.radius, self.mass, self.meta['kernels'], self.spkid,
                   self.error_ra, self.error_dec)
               )
        if hasattr(self, 'offset'):
            out += 'Offset applied: RA*cosDEC: {:.4f}; DEC: {:.4f}\n'.format(
                self.offset.d_lon_coslat.to(u.arcsec),
                self.offset.d_lat.to(u.arcsec))
        return out
