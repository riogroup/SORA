"""Base ephemeris support shared by SORA ephemeris classes."""

import warnings

import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.time import Time

from sora.config import input_tests
from sora.config.decorators import deprecated_function


class BaseEphem:
    """Base class with shared metadata and helper methods for ephemerides.

    Parameters
    ----------
    name : `str`, optional
        Name of the object.

    spkid : `str`, `int`, optional
        SPK identifier of the object.

    error_ra : `int`, `float`, `astropy.units.Quantity`, optional
        Ephemeris RA*cosDEC uncertainty. Values without units are interpreted
        in arcsec.

    error_dec : `int`, `float`, `astropy.units.Quantity`, optional
        Ephemeris DEC uncertainty. Values without units are interpreted in
        arcsec.

    radius : `int`, `float`, `astropy.units.Quantity`, optional
        Object radius. Values without units are interpreted in km.

    mass : `int`, `float`, optional, default=0
        Object mass, in kg.

    H : `int`, `float`, optional
        Object absolute magnitude.

    G : `int`, `float`, optional
        Object phase slope.
    """

    def __init__(self, name=None, spkid=None, **kwargs):
        """Initialize shared ephemeris attributes.

        Parameters
        ----------
        name : `str`, optional
            Name of the object.

        spkid : `str`, `int`, optional
            SPK identifier of the object.

        **kwargs
            Optional ephemeris uncertainty and physical parameters accepted by
            `BaseEphem`.
        """
        # remove 'H', 'G', 'mass' and 'radius' from allowed kwargs and docstring for v1.0
        input_tests.check_kwargs(kwargs, allowed_kwargs=['error_dec', 'error_ra', 'H', 'G', 'mass', 'radius'])
        #
        self._shared_with = {'body': {}, 'occultation': {}}
        self.name = name
        if spkid:
            self.spkid = spkid
            self.code = self.spkid  # remove this line for v1.0
        self.error_ra = u.Quantity(kwargs.get('error_ra', 0), unit=u.arcsec)
        if self.error_ra < 0:
            warnings.warn("Error in RA cannot be negative. Using absolute value.")
            self.error_ra = np.absolute(self.error_ra)
        self.error_dec = u.Quantity(kwargs.get('error_dec', 0), unit=u.arcsec)
        if self.error_dec < 0:
            warnings.warn("Error in DEC cannot be negative. Using absolute value.")
            self.error_dec = np.absolute(self.error_dec)
        self.offset = (0, 0)
        # start of block removal for v1.0
        if 'radius' in kwargs:
            self.radius = kwargs['radius'] * u.km
        self.mass = kwargs.get('mass', 0.0) * u.kg
        if 'H' in kwargs:
            self.H = kwargs['H']
        if 'G' in kwargs:
            self.G = kwargs['G']
        # end of block removal for v1.0

    @property
    def spkid(self):
        """`str` : SPK identifier of the object."""
        if 'spkid' in self._shared_with['body']:
            return self._shared_with['body']['spkid']
        elif hasattr(self, '_spkid'):
            return self._spkid
        else:
            raise AttributeError('{} does not have spkid'.format(self.__class__.__name__))

    @spkid.setter
    def spkid(self, value):
        """Set the object's SPK identifier."""
        if 'spkid' in self._shared_with['body']:
            raise AttributeError('When {} is associated to a Body object, spkid must be given to the Body'
                                 ' object.'.format(self.__class__.__name__))
        self._spkid = str(int(value))

    # Start of block removal for v1.0
    @property
    def radius(self):
        """`astropy.units.Quantity` : Object radius, in km."""
        if 'radius' in self._shared_with['body']:
            return self._shared_with['body']['radius']
        elif hasattr(self, '_radius'):
            return self._radius
        else:
            raise AttributeError('{} does not have radius'.format(self.__class__.__name__))

    @radius.setter
    def radius(self, value):
        """Set the object radius."""
        if 'radius' in self._shared_with['body']:
            raise AttributeError('When {} is associated to a Body object, radius must be given to the Body'
                                 ' object.'.format(self.__class__.__name__))
        self._radius = u.Quantity(value, unit=u.km)

    @property
    def H(self):
        """`float` : Object absolute magnitude."""
        if 'H' in self._shared_with['body']:
            return self._shared_with['body']['H'].value
        elif hasattr(self, '_H'):
            return self._H
        else:
            raise AttributeError('{} does not have H'.format(self.__class__.__name__))

    @H.setter
    def H(self, value):
        """Set the object absolute magnitude."""
        if 'H' in self._shared_with['body']:
            raise AttributeError('When {} is associated to a Body object, H must be given to the Body'
                                 ' object.'.format(self.__class__.__name__))
        self._H = float(value)

    @property
    def G(self):
        """`float` : Object phase slope."""
        if 'G' in self._shared_with['body']:
            return self._shared_with['body']['G']
        elif hasattr(self, '_G'):
            return self._G
        else:
            raise AttributeError('{} does not have G'.format(self.__class__.__name__))

    @G.setter
    def G(self, value):
        """Set the object phase slope."""
        if 'G' in self._shared_with['body']:
            raise AttributeError('When {} is associated to a Body object, G must be given to the Body'
                                 ' object.'.format(self.__class__.__name__))
        self._G = float(value)

    def apparent_magnitude(self, time):
        """Calculates the object's apparent magnitude.

        Parameters
        ----------
        time : `str`, `astropy.time.Time`
            Reference time to calculate the object's apparent magnitude.
            It can be a string in the ISO format (yyyy-mm-dd hh:mm:ss.s) or an astropy Time object.

        Returns
        -------
        ap_mag : `float`, `list`
            Object apparent magnitude. A list is returned when multiple times
            are requested and the value is obtained from Horizons.
        """
        from sora.body.utils import apparent_magnitude
        from astroquery.jplhorizons import Horizons
        from astropy.coordinates import get_sun

        time = Time(time)

        if getattr(self, 'H', None) is None or getattr(self, 'G', None) is None:
            search_name = self._shared_with['body'].get('search_name', self.name)
            id_type = getattr(self, 'id_type', 'majorbody')
            id_type = self._shared_with['body'].get('id_type', id_type)
            obj = Horizons(id=search_name, id_type=id_type, location='geo', epochs=time.jd)
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
            obs_sun = get_sun(time)
            sun_obj = SkyCoord(obs_obj.cartesian - obs_sun.cartesian)
            sun_obj.representation_type = 'spherical'

            # Calculates the phase angle between the 2-vectors
            unit_vector_1 = -obs_obj.cartesian.xyz / np.linalg.norm(obs_obj.cartesian.xyz)
            unit_vector_2 = -sun_obj.cartesian.xyz / np.linalg.norm(sun_obj.cartesian.xyz)
            dot_product = np.dot(unit_vector_1, unit_vector_2)
            phase = np.arccos(dot_product).to(u.deg).value

            return apparent_magnitude(self.H, self.G, obs_obj.distance.to(u.AU).value,
                                      sun_obj.distance.to(u.AU).value, phase)

    @deprecated_function(message="Please use get_pole_position_angle from Body object")
    def get_pole_position_angle(self, pole, time):
        """Returns the pole position and aperture angles.

        Parameters
        ----------
        pole : `str`, `astropy.coordinates.SkyCoord`
            Coordinate of the object pole ICRS.

        time : `str`, `astropy.time.Time`
            Time from which to calculate the position. It can be a string
            in the ISO format (yyyy-mm-dd hh:mm:ss.s) or an astropy Time object.

        Returns
        -------
        position_angle : `float`
            Position angle of the object pole, in degrees.

        aperture_angle : `float`
            Aperture angle of the object pole, in degrees.
        """
        time = Time(time)
        if type(pole) == str:
            pole = SkyCoord(pole, unit=(u.hourangle, u.deg))
        obj = self.get_position(time)
        position_angle = obj.position_angle(pole).value * u.rad
        aperture_angle = np.arcsin(
            -(np.sin(pole.dec) * np.sin(obj.dec) +
              np.cos(pole.dec) * np.cos(obj.dec) * np.cos(pole.ra - obj.ra))
        )
        return position_angle.to('deg'), aperture_angle.to('deg')

    # End of block removal for v1.0

    def get_ksi_eta(self, time, star):
        """Returns the object's projected position relative to a star.

        Returns the projected position (orthographic projection) of the object
        in the tangent sky plane relative to a star.

        Parameters
        ----------
        time : `str`, `astropy.time.Time`
            Reference time to calculate the object position. It can be a string
            in the ISO format (yyyy-mm-dd hh:mm:ss.s) or an astropy Time object.

        star : `str`, `astropy.coordinates.SkyCoord`
            Coordinate of the star in the same reference frame as the ephemeris.

        Returns
        -------
        ksi, eta : `float`
            Projected position (orthographic projection) of the object in the
            tangent sky plane relative to a star.
            ``ksi`` is in the East-West direction (East positive).
            ``eta`` is in the North-South direction (North positive).
        """
        from astropy.coordinates import SkyOffsetFrame
        time = Time(time)
        if type(star) == str:
            star = SkyCoord(star, unit=(u.hourangle, u.deg))
        coord = self.get_position(time)
        target = coord.transform_to(SkyOffsetFrame(origin=star))
        da = target.cartesian.y
        dd = target.cartesian.z
        return da.to(u.km).value, dd.to(u.km).value

    def add_offset(self, da_cosdec, ddec):
        """Sets the offset applied to the ephemeris.

        Parameters
        ----------
        da_cosdec : `int`, `float`
            Delta_alpha_cos_delta, in mas.

        ddec : `int`, `float`
            Delta_delta, in mas.
        """
        self.offset = (da_cosdec, ddec)

    @property
    def offset(self):
        """`astropy.coordinates.SphericalCosLatDifferential` : Ephemeris offset."""
        return self._offset

    @offset.setter
    def offset(self, value):
        """Set the ephemeris offset from RA*cosDEC and DEC values."""
        from astropy.coordinates import SphericalCosLatDifferential
        dadc, dd = value
        self._offset = SphericalCosLatDifferential(dadc * u.mas, dd * u.mas, 0.0 * u.km)

    def to_file(self, time, namefile=None):
        """Saves the ephemerides to a file.

        Ephemeris will be saved starting one hour before the central time
        until one hour after it, with a step of one minute.

        Note
        ----
            This file can be used as an input for ``EphemPlanete()``.

        Parameters
        ----------
        time : `str`, `astropy.time.Time`
            Central time to be saved.

        namefile : `str`
            Filename to save. If not given, a default filename is generated
            from the object name.
        """
        if namefile is None:
            namefile = 'Ephem_' + self.name.replace(' ', '_') + '.dat'
        time = input_tests.test_attr(time, Time, 'time')
        time_output = time + np.arange(-60, 61, 1) * u.min
        ephem_output = self.get_position(time_output.utc)
        array_output = np.array([time_output.utc.jd, ephem_output.ra.deg, ephem_output.dec.deg,
                                 ephem_output.distance.au]).T
        np.savetxt(namefile, array_output, delimiter='    ', fmt='%.14f')

    def __str__(self):
        """String representation of the Ephem Class.
        """
        out = ("----------- Ephemeris -----------\n"
               "\n{}: {{ephem_info}} (SPKID={})\n"
               "Ephem Error: RA*cosDEC: {:.3f}; DEC: {:.3f}\n".format(
                   self.__class__.__name__, getattr(self, 'spkid', ''), self.error_ra, self.error_dec)
               )
        if hasattr(self, 'offset'):
            out += 'Offset applied: RA*cosDEC: {:.4f}; DEC: {:.4f}\n'.format(
                self.offset.d_lon_coslat.to(u.arcsec),
                self.offset.d_lat.to(u.arcsec))
        return out
