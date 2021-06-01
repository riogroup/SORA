import astropy.units as u
from astropy.coordinates import SkyCoord, EarthLocation, GCRS, AltAz
from astropy.time import Time

from sora.config import input_tests
from .utils import search_code_mpc

__all__ = ['Observer', 'Spacecraft']


class Observer:
    """Defines the observer object.

    Attributes
    ----------
    name : `str`
        Name for the Observer. Observer is uniquely defined (name must be
        different for each observer).

    code : `str`
        The IAU code for SORA to search for its coordinates in MPC database.

    site : `astropy.coordinates.EarthLocation`
        User provides an EarthLocation object.

    lon : `str`, `float`
        The Longitude of the site in degrees.
        Positive to East. Range (0 to 360) or (-180 to +180).
        User can provide in degrees (`float`) or hexadecimal (`string`).

    lat : `str`, `float`
        The Latitude of the site in degrees.
        Positive North. Range (+90 to -90).
        User can provide in degrees (float) or hexadecimal (string).

    height : `int`, `float`
        The height of the site in meters above see level.

    ephem : `str`, `list`
        The ephemeris used to locate the observer on space.
        It can be "horizons" to use horizons or a list of kernels


    Examples
    --------
    User can provide one of the following to define an observer:

    - If user will use the MPC name for the site:

    >>> Observer(code)

    - If user wants to use a different name from the MPC database:

    >>> Observer(name, code)

    - If user wants to use an EarthLocation value:

    >>> from astropy.coordinates import EarthLocation
    >>> EarthLocation(lon, lat, height)
    >>> Observer(name, site)

    - If user wants to give site coordinates directly:

    >>> Observer(name, lon, lat, height)

    """

    def __init__(self, **kwargs):

        input_tests.check_kwargs(kwargs, allowed_kwargs=['code', 'height', 'lat', 'lon', 'name', 'site', 'ephem'])
        self.__name = kwargs.get('name', '')
        if 'code' in kwargs and any(i in kwargs for i in ['lon', 'lat', 'height']):
            raise ValueError("The Observer object is instantiated with IAU code or coordinates, not both.")
        if 'code' in kwargs:
            self.code = kwargs['code']
            try:
                name, self.site = search_code_mpc()[self.code]
                self.__name = kwargs.get('name', name)
            except:
                raise ValueError('code {} could not be located in MPC database'.format(self.code))
        elif 'site' in kwargs:
            self.site = input_tests.test_attr(kwargs['site'], EarthLocation, 'site')
        elif all(i in kwargs for i in ['lon', 'lat']):
            self.site = EarthLocation(kwargs['lon'], kwargs['lat'], kwargs.get('height', 0.0))
        else:
            raise ValueError('Input parameters could not be determined')
        self.ephem = kwargs.get('ephem', 'horizons')

    def get_ksi_eta(self, time, star):
        """Calculates relative position to star in the orthographic projection.

        Parameters
        ----------
        time : `str`, `astropy.time.Time`
            Reference time to calculate the position. It can be a string in the
            format ``'yyyy-mm-dd hh:mm:ss.s'`` or an astropy `Time` object.

        star : `str`, `astropy.coordinates.SkyCoord`
            The coordinate of the star in the same reference frame as the
            ephemeris. It can be a string in the format ``'hh mm ss.s +dd mm ss.ss'``
            or an astropy `SkyCoord` object.


        Returns
        -------
        ksi, eta : `float`
            On-sky orthographic projection of the observer relative to a star.
            ``ksi`` is in the North-South direction (North positive).
            ``eta`` is in the East-West direction (East positive).
        """
        from astropy.coordinates.matrix_utilities import rotation_matrix

        time = input_tests.test_attr(time, Time, 'time')
        try:
            star = SkyCoord(star, unit=(u.hourangle, u.deg))
        except:
            raise ValueError('star is not an astropy object or a string in the format "hh mm ss.s +dd mm ss.ss"')

        itrs = self.site.get_itrs(obstime=time)
        gcrs = itrs.transform_to(GCRS(obstime=time))
        rz = rotation_matrix(star.ra, 'z')
        ry = rotation_matrix(-star.dec, 'y')

        cp = gcrs.cartesian.transform(rz).transform(ry)
        return cp.y.to(u.km).value, cp.z.to(u.km).value

    def sidereal_time(self, time, mode='local'):
        """Calculates the Apparent Sidereal Time at a reference time.

        Parameters
        ----------
        time : `str`, `astropy.time.Time`
            Reference time to calculate sidereal time.It can be a string
            in the ISO format (yyyy-mm-dd hh:mm:ss.s) or an astropy Time object.

        mode : `str`
            Local or greenwich time.
            If mode set ``'local'`` calculates the sidereal time for the
            coordinates of this object.
            If mode set ``'greenwich'`` calculates the Greenwich Apparent
            Sidereal Time.


        Returns
        -------
        sidereal_time
            An Astropy Longitude object with the Sidereal Time.
        """
        # return local or greenwich sidereal time
        time = input_tests.test_attr(time, Time, 'time')
        time.location = self.site
        if mode == 'local':
            return time.sidereal_time('apparent')
        elif mode == 'greenwich':
            return time.sidereal_time('apparent', 'greenwich')
        else:
            raise ValueError('mode must be "local" or "greenwich"')

    def altaz(self, time, coord):
        """Calculates the Altitude and Azimuth at a reference time for a coordinate.

        Parameters
        ----------
        time : `str`, `astropy.time.Time`
            Reference time to calculate the sidereal time. It can be a string
            in the ISO format (yyyy-mm-dd hh:mm:ss.s) or an astropy Time object.

        coord : `str`, `astropy.coordinates.SkyCoord`
            Coordinate of the target ICRS.


        Returns
        -------
        altitude : `float`
            Object altitude in degrees.

        azimuth : `float`
            Object azimuth in degrees.
        """
        time = input_tests.test_attr(time, Time, 'time')
        if type(coord) == str:
            coord = SkyCoord(coord, unit=(u.hourangle, u.deg))
        ephem_altaz = coord.transform_to(AltAz(obstime=time, location=self.site))
        return ephem_altaz.alt.deg, ephem_altaz.az.deg

    def get_vector(self, time, origin='barycenter'):
        """Return the vector Origin -> Observer in the ICRS

        Parameters
        ----------
        time : `str`, `astropy.time.Time`
            Reference time to calculate the object position. It can be a string
            in the ISO format (yyyy-mm-dd hh:mm:ss.s) or an astropy Time object.

        origin : `str`
            Origin of vector. It can be 'barycenter' or 'geocenter'.

        Returns
        -------
        coord : `astropy.coordinates.SkyCoord`
            Astropy SkyCoord object with the vector origin -> observer at the given time.
        """
        from sora.ephem.utils import ephem_horizons, ephem_kernel

        itrs = self.site.get_itrs(obstime=time)
        gcrs = itrs.transform_to(GCRS(obstime=time))

        if self.ephem == 'horizons':
            vector = ephem_horizons(time=time, target=self.spkid, observer=origin, id_type='majorbody', output='vector')
        else:
            vector = ephem_kernel(time=time, target=self.spkid, observer=origin, kernels=self.ephem, output='vector')

        topo = SkyCoord(vector.cartesian + gcrs.cartesian, representation_type='cartesian')
        if not topo.isscalar and len(topo) == 1:
            topo = topo[0]
        return topo

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

    @property
    def spkid(self):
        return '399'

    def __repr__(self):
        """String representation of the Observer Class
        """
        return '<{}: {}>'.format(self.__class__.__name__, self.name)

    def __str__(self):
        """ String representation of the Observer class
        """
        out = ('Site: {}\n'
               'Geodetic coordinates: Lon: {}, Lat: {}, height: {:.3f}'.format(
                   self.name, self.site.lon.__str__(), self.site.lat.__str__(),
                   self.site.height.to(u.km))
               )
        return out


class Spacecraft:
    """Defines a spacecraft observer object.

    Attributes
    ----------
    name : `str`
        Name for the Observer. Observer is uniquely defined (name must be
        different for each observer).

    spkid : `str`, required
        `spkid` of the targeting object.

    ephem : `str`, `list`
        The ephemeris used to locate the observer on space.
        It can be "horizons" to use horizons or a list of kernels

    """

    def __init__(self, name, spkid, ephem='horizons'):
        self._name = name
        self._spkid = spkid
        self._ephem = ephem

    @property
    def name(self):
        return self._name

    @property
    def spkid(self):
        return self._spkid

    @property
    def ephem(self):
        return self._ephem

    def get_vector(self, time, origin='barycenter'):
        """Return the vector Origin -> Observer in the ICRS

        Parameters
        ----------
        time : `str`, `astropy.time.Time`
            Reference time to calculate the object position. It can be a string
            in the ISO format (yyyy-mm-dd hh:mm:ss.s) or an astropy Time object.

        origin : `str`
            Origin of vector. It can be 'barycenter' or 'geocenter'.

        Returns
        -------
        coord : `astropy.coordinates.SkyCoord`
            Astropy SkyCoord object with the vector origin -> observer at the given time.
        """
        from sora.ephem.utils import ephem_horizons, ephem_kernel
        time = Time(time)

        if self.ephem == 'horizons':
            vector = ephem_horizons(time=time, target=self.spkid, observer=origin, id_type='majorbody', output='vector')
        else:
            vector = ephem_kernel(time=time, target=self.spkid, observer=origin, kernels=self.ephem, output='vector')
        return vector

    def __repr__(self):
        """String representation of the Spacecraft Class
        """
        return '<{}: {}>'.format(self.__class__.__name__, self.name)

    def __str__(self):
        """ String representation of the Spacecraft class
        """
        out = ['Spacecraft: {} (spkid={})'.format(self.name, self.spkid),
               'Positions from {}'.format(str(self.ephem))]
        return '\n'.join(out)
