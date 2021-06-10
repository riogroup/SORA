import warnings

import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord, SkyOffsetFrame, ICRS
from astropy.time import Time

from sora.config.decorators import deprecated_alias
from .meta import BaseEphem

__all__ = ['EphemHorizons', 'EphemJPL', 'EphemKernel', 'EphemPlanete']


class EphemPlanete(BaseEphem):
    """Class used to simulate former Fortran programs `ephem_planete` and
    `fit_d2_ksi_eta`.

    Attributes
    ----------
    ephem : `file`, required
        Input file with JD (UTC), geocentric RA (deg), DEC (deg), and
        distance (AU).

    name : `str`, optional, default=None
        Name of the object to search in the JPL database.

    radius : `int`, `float`, optional, default: online database
        Object radius, in km.

    error_ra : `int`, `float`, optional, default: online database
        Ephemeris RA*cosDEC error, in arcsec.

    error_dec : `int`, `float`, optional, default: online database
        Ephemeris DEC error, in arcsec.

    mass : `int`, `float`, optional. default=0
        Object mass, in kg.

    H : `int`, `float`, optional, default=NaN
        Object absolute magnitude.

    G : `int`, `float`, optional, default=NaN
        Object phase slope.

    """

    def __init__(self, ephem, name=None, spkid=None, **kwargs):

        super().__init__(name=name, spkid=spkid, **kwargs)
        data = np.loadtxt(ephem, unpack=True)
        self.time = Time(data[0], format='jd')
        self.ephem = SkyCoord(data[1] * u.deg, data[2] * u.deg, data[3] * u.AU)
        self.__reftime = self.time[0]
        self.min_time = Time(data[0].min(), format='jd')
        self.max_time = Time(data[0].max(), format='jd')

    def get_position(self, time, observer='geocenter'):
        """Returns the object's geocentric position.

        Parameters
        ----------
        time : `str`, `astropy.time.Time`
            Reference time to calculate the object position. It can be a string
            in the ISO format (yyyy-mm-dd hh:mm:ss.s) or an astropy Time object.

        observer : `any`
            This parameter is present in EphemPlanete for compatibility with the
            remaining ephem classes. The returned positions are based on the given
            ephemeris despite the observer.

        Returns
        -------
        coord : `astropy.coordinates.SkyCoord`
            Astropy SkyCoord object with the object coordinates at the given time.
        """
        ksi, eta = self.get_ksi_eta(time=time) * u.km
        distance = self.ephem.distance.mean()
        off_ra = np.arctan2(ksi, distance)
        off_dec = np.arctan2(eta, distance)
        coord_frame = SkyOffsetFrame(origin=self.star)
        pos = SkyCoord(lon=off_ra, lat=off_dec, distance=distance, frame=coord_frame)
        return pos.icrs

    @deprecated_alias(log='verbose')  # remove this line in v1.0
    def fit_d2_ksi_eta(self, star, verbose=True):
        """Fits the projected position (orthographic projection) of the object in
        the tangent sky plane relative to a star.

        Parameters
        ----------
        star : `str`, `astropy.coordinates.SkyCoord`
            The coordinate of the star in the same reference frame as the ephemeris.

        verbose : `bool`, optional, default=True
            Enable log printing.
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
        if verbose:
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
        """Returns the projected position (orthographic projection) of the object
        in the tangent sky plane relative to a star.

        Parameters
        ----------
        time : `str`, `astropy.time.Time`, required
            Reference time to calculate the position. It can be a string
            in the ISO format (yyyy-mm-dd hh:mm:ss.s) or an astropy Time object.

        star : `str`, `astropy.coordinates.SkyCoord`, optional, default=None
            The coordinate of the star in the same reference frame as the ephemeris.

        Returns
        -------
        ksi, eta : `float` array
            Projected position (orthographic projection) of the object in the
            tangent sky plane relative to a star.
            ``ksi`` is in the North-South direction (North positive).
            ``eta`` is in the East-West direction (East positive).
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
    """Obtains the ephemeris from Horizons/JPL service.

    Note
    ----
    Web tool URL: https://ssd.jpl.nasa.gov/horizons.cgi


    Attributes
    ----------
    name : `str`, required
        Name of the object to search in the JPL database.

    id_type: `str`, default='smallbody'
        Type of object options: ``smallbody``, ``majorbody`` (planets but
        also anything that is not a small body), ``designation``, ``name``,
        ``asteroid_name``, ``comet_name``, ``id`` (Horizons id number), or
        ``smallbody`` (find the closest match under any id_type).

    radius : `int`, `float`, default: online database
        Object radius, in km.

    error_ra : `int`, `float`, default: online database
        Ephemeris RA*cosDEC error, in arcsec.

    error_dec : `int`, `float`, default: online database
        Ephemeris DEC error, in arcsec.

    mass : `int`, `float`, default=0
        Object mass, in kg.

    H : `int`, `float`, default=NaN
        Object absolute magnitude.

    G : `int`, `float`, default=NaN
        Object phase slope.

    """

    def __init__(self, name, id_type='smallbody', spkid=None, **kwargs):

        super().__init__(name=name, spkid=spkid, **kwargs)
        self.id_type = id_type
        _ = self.get_position(Time.now())  # test if Horizons can proceed for this object

    def get_position(self, time, observer='geocenter'):
        """Returns the ICRS position of the object for observer.

        Parameters
        ----------
        time : `str`, `astropy.time.Time`
            Reference time to calculate the object position. It can be a string
            in the ISO format (yyyy-mm-dd hh:mm:ss.s) or an astropy Time object.

        observer : `str`, `sora.Observer`, `sora.Spacecraft`
            IAU code of the observer (must be present in given list of kernels),
            a SORA observer object or a string: ['geocenter', 'barycenter']

        Returns
        -------
        coord : `astropy.coordinates.SkyCoord`
            Astropy SkyCoord object with the object coordinates at the given time.
        """
        from .utils import ephem_horizons

        coord = ephem_horizons(time=time, target=self.name, observer=observer, id_type=self.id_type, output='ephemeris')
        if hasattr(self, 'offset'):
            pos_frame = SkyOffsetFrame(origin=coord)
            new_pos = SkyCoord(lon=self.offset.d_lon_coslat, lat=self.offset.d_lat,
                               distance=coord.distance, frame=pos_frame)
            coord = new_pos.transform_to(ICRS)
        if not coord.isscalar and len(coord) == 1:
            coord = coord[0]
        return coord

    def __str__(self):
        """ String representation of the EphemHorizons Class.
        """
        out = super().__str__().format(ephem_info='Ephemeris are downloaded from Horizons website')
        return out


# remove this block for v1.0
class EphemJPL(EphemHorizons):
    def __init__(self, name, id_type='smallbody', spkid=None, **kwargs):
        warnings.warn('EphemJPL is deprecated and will be removed in v1.0. Please use EphemHorizons.')
        super().__init__(name=name, id_type=id_type, spkid=spkid, **kwargs)

    __init__.__doc__ = EphemHorizons.__init__.__doc__
# end of block removal for v1.0


class EphemKernel(BaseEphem):
    """Gets the ephemeris from BSP kernels.

    Parameters
    ----------
    name : `str`,  optional, default=None
        Name of the object to search in the JPL database.

    spkid : `str`, required
        `spkid` of the targeting object. Former 'code' (v0.1).

    kernels : `list`, required
        List of paths for kernels files.

    radius : `int`, `float`, optional, default: online database
        Object radius, in km.

    error_ra : `int`, `float`, optional, default: online database
        Ephemeris RA*cosDEC error, in arcsec .

    error_dec : `int`, `float`, optional, default: online database
        Ephemeris DEC error, in arcsec.

    mass : `int`, `float`, optional, default=0
        Object Mass, in kg.

    H : `int`, `float`, optional, default=NaN
        Object Absolute Magnitude.

    G : `int`, `float`, optional, default=NaN
        Object Phase slope.

    """

    @deprecated_alias(code='spkid')  # remove this line for v1.0
    def __init__(self, kernels, spkid, name=None, **kwargs):
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

    def get_position(self, time, observer='geocenter'):
        """Returns the object geocentric position.

        Parameters
        ----------
        time : `str`, `astropy.time.Time`
            Reference time to calculate the object position. It can be a string
            in the ISO format (yyyy-mm-dd hh:mm:ss.s) or an astropy Time object.

        observer : `str`, `sora.Observer`, `sora.Spacecraft`
            IAU code of the observer (must be present in given list of kernels),
            a SORA observer object or a string: ['geocenter', 'barycenter']

        Returns
        -------
        coord : `astropy.coordinates.SkyCoord`
            Astropy SkyCoord object with the object coordinates at the given time.
        """
        from .utils import ephem_kernel
        pos = ephem_kernel(time=time, target=self.spkid, observer=observer, kernels=self.__kernels)
        if hasattr(self, 'offset'):
            pos_frame = SkyOffsetFrame(origin=pos)
            new_pos = SkyCoord(lon=self.offset.d_lon_coslat, lat=self.offset.d_lat,
                               distance=pos.distance, frame=pos_frame)
            pos = new_pos.transform_to(ICRS)
        pos = SkyCoord(ra=pos.ra, dec=pos.dec, distance=pos.distance)
        if not pos.isscalar and len(pos) == 1:
            pos = pos[0]
        return pos

    def __str__(self):
        """ String representation of the EphemKernel Class.
        """
        out = super().__str__().format(ephem_info=self.meta['kernels'])
        return out
