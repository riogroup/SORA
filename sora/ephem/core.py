import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord, SkyOffsetFrame, ICRS
from astropy.time import Time

from sora.config.decorators import deprecated_alias
from .meta import BaseEphem

__all__ = ['EphemHorizons', 'EphemJPL', 'EphemKernel', 'EphemPlanete']


class EphemPlanete(BaseEphem):
    def __init__(self, ephem, name=None, spkid=None, **kwargs):
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
        super().__init__(name=name, spkid=spkid, **kwargs)
        data = np.loadtxt(ephem, unpack=True)
        self.time = Time(data[0], format='jd')
        self.ephem = SkyCoord(data[1] * u.deg, data[2] * u.deg, data[3] * u.AU)
        self.__reftime = self.time[0]
        self.min_time = Time(data[0].min(), format='jd')
        self.max_time = Time(data[0].max(), format='jd')

    def get_position(self, time):
        """ Returns the geocentric position of the object.

        Parameters:
            time (str, Time):Time from which to calculate the position.

        Returns:
            coord (SkyCoord): Astropy SkyCoord object with the coordinate at given time
        """
        ksi, eta = self.get_ksi_eta(time=time) * u.km
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
        if hasattr(self, 'star') and self.star.to_string('hmsdms', precision=5) == star.to_string('hmsdms',
                                                                                                  precision=5):
            return
        self.star = star
        target = self.ephem.transform_to(SkyOffsetFrame(origin=star))
        da = target.cartesian.y
        dd = target.cartesian.z
        dt = (self.time - self.__reftime) / (self.max_time - self.min_time)

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
            k = ksi((time - self.__reftime) / (self.max_time - self.min_time))
            e = eta((time - self.__reftime) / (self.max_time - self.min_time))
            if hasattr(self, 'offset'):
                dist = self.ephem[0].distance.to(u.km).value
                dksi = np.sin(self.offset.d_lon_coslat) * dist
                deta = np.sin(self.offset.d_lat) * dist
                k = k + dksi
                e = e + deta
            return k.value, e.value
        else:
            raise ValueError('A "star" parameter is missing. Please run fit_d2_ksi_eta first.')

    def __str__(self):
        """ String representation of the EphemPlanete Class.
        """
        validity = 'Valid from {} until {}'.format(self.min_time.iso, self.max_time.iso)
        out = super().__str__().format(ephem_info=validity)
        if hasattr(self, 'star'):
            out += ("\nFitted ephemeris position relative to star coordinate {}\n"
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
            dt = (self.time - self.__reftime) / (self.max_time - self.min_time)
            ksi = np.poly1d(self.ksi)
            eta = np.poly1d(self.eta)
            dksi = da.to(u.km).value - ksi(dt)
            deta = dd.to(u.km).value - eta(dt)
            rmsk = np.sqrt(np.mean(np.square(dksi)))
            rmse = np.sqrt(np.mean(np.square(deta)))
            out += 'Residual RMS: ksi={:.3f} km, eta={:.3f} km\n'.format(rmsk, rmse)
        return out


class EphemHorizons(BaseEphem):
    def __init__(self, name, id_type='smallbody', spkid=None, **kwargs):
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
        super().__init__(name=name, spkid=spkid, **kwargs)
        self.id_type = id_type

    def get_position(self, time):
        """ Returns the geocentric position of the object.

        Parameters:
            time (str, Time): Reference time to calculate the position.

        Returns:
            coord (SkyCoord): Astropy SkyCoord object with the object coordinates at the given time
        """
        from astroquery.jplhorizons import Horizons
        from astropy.table import vstack

        time = Time(time)
        if not time.isscalar:
            s = len(time)

            def calc_ephem(i):
                n = 50 * i
                k = n + 50
                if k > s:
                    k = s
                ob = Horizons(id=self.name, id_type=self.id_type, location='geo', epochs=time[n:k].jd)
                return ob.ephemerides(extra_precision=True)

            plus = 0
            if s % 50 > 0:
                plus = 1
            eph = vstack([calc_ephem(j) for j in range(s // 50 + 1 * plus)])
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

    def __str__(self):
        """ String representation of the EphemHorizons Class.
        """
        out = super().__str__().format(ephem_info='Ephemeris are downloaded from Horizons website')
        return out


# remove this block for v1.0
class EphemJPL(EphemHorizons):
    pass
# end of block removal for v1.0


class EphemKernel(BaseEphem):
    @deprecated_alias(code='spkid')  # remove this line for v1.0
    def __init__(self, kernels, spkid, name=None, **kwargs):
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
        import spiceypy as spice

        super().__init__(name=name, spkid=spkid, **kwargs)
        self.meta = {}
        kerns = []
        for arg in kernels:
            spice.furnsh(arg)
            kerns.append(arg.split('/')[-1].split('.')[0].upper())
        spice.kclear()
        self.meta['kernels'] = '/'.join(kerns)
        self.__kernels = kernels

    def get_position(self, time):
        """ Returns the object geocentric position.

        Parameters:
            time (str, Time): Reference time to calculate the object position.

        Returns:
            coord (SkyCoord): Astropy SkyCoord object with the object coordinates at the given time
        """
        from .utils import ephem_kernel
        pos = ephem_kernel(time, self.spkid, '399', self.__kernels)
        if hasattr(self, 'offset'):
            pos_frame = SkyOffsetFrame(origin=pos)
            new_pos = SkyCoord(lon=self.offset.d_lon_coslat, lat=self.offset.d_lat,
                               distance=pos.distance, frame=pos_frame)
            pos = new_pos.transform_to(ICRS)
        return pos

    def __str__(self):
        """ String representation of the EphemKernel Class.
        """
        out = super().__str__().format(ephem_info=self.meta['kernels'])
        return out
