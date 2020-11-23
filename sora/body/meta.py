from sora.ephem import EphemPlanete, EphemKernel, EphemJPL, EphemHorizons
import astropy.units as u
import astropy.constants as const
from astropy.coordinates import SkyCoord, Longitude, Latitude
import numpy as np
import warnings


__all__ = ['PhysicalData']


class PhysicalData(u.Quantity):
    """ Define PhysicalData with uncertainty, reference and notes.
    It inherits from astropy.units.quantity.Quantity()

    Parameters:
        name (str): The name representing the corresponding physical parameter.
        value (number, `~numpy.ndarray`, `Quantity` object (sequence), str):
            The numerical value of this quantity in the units given by unit.  If a
            `Quantity`  (or any other valid object with a ``unit`` attribute),
            creates a new `Quantity` object, converting to `unit` units as needed.
            If a string, it is converted to a number or `Quantity`, depending on
            whether a unit is present.
        uncertainty (number, `~numpy.ndarray`, `Quantity` object (sequence), str):
            The numerical value of this quantity in the units given by unit.  If a
            `Quantity` (or any other valid object with a ``unit`` attribute),
            creates a new `Quantity` object, converting to `unit` units as needed.
            If a string, it is converted to a number or `Quantity`, depending on
            whether a unit is present. Default = 0.0
        reference (str): A string stating the reference for the parameter value.
            Default = "User"
        notes (str): Any other important information about the physical parameter.
            Default = ""
        unit: `~astropy.units.UnitBase` instance, str
            An object that represents the unit associated with the input value.
            Must be an `~astropy.units.UnitBase` object or a string parseable by
            the :mod:`~astropy.units` package. Default = "dimensionless"
        raise_error (bool): If value=None, the function raise an error if
            raise_error=True, else value is redefined to NaN.
    """

    def __new__(cls, name, value, uncertainty=0.0, reference="User", notes="", unit=u.dimensionless_unscaled,
                raise_error=False):
        given_unit = unit
        if value is None:
            if raise_error:
                raise TypeError("The value must be a valid PhysicalData type. Given: {}".format(value))
            else:
                value = np.nan
        elif isinstance(value, u.core.CompositeUnit):
            unit = u.dimensionless_unscaled
            for i in np.arange(len(value.bases)):
                unit = unit*(value.bases[i]**value.powers[i])
            value = value.scale
        elif isinstance(value, (PhysicalData, u.Quantity)):
            unit = value.unit
            value = value.value
        if not np.isscalar(value):
            print(value)
            raise ValueError('Given value must be a scalar. Given: {}'.format(value))
        if not unit.is_equivalent(given_unit):
            raise ValueError('{} is not equivalent to {}'.format(unit, given_unit))
        physdata = super().__new__(cls, value=value, unit=unit)
        physdata = physdata.to(given_unit)
        physdata.name = name
        physdata.uncertainty = uncertainty
        physdata.reference = reference
        physdata.notes = notes
        return physdata

    @property
    def uncertainty(self):
        return self._uncertainty

    @uncertainty.setter
    def uncertainty(self, value):
        unit = self.unit
        given_unit = unit
        if isinstance(value, u.core.CompositeUnit):
            unit = u.dimensionless_unscaled
            for i in np.arange(len(value.bases)):
                unit = unit*(value.bases[i]**value.powers[i])
            value = value.scale
        elif value is None:
            value = 0.0
        elif isinstance(value, u.Quantity):
            unit = value.unit
            value = value.value
        if not np.isscalar(value):
            print(value)
            raise ValueError('Given value must be a scalar. Given: {}'.format(value))
        if not unit.is_equivalent(self.unit):
            raise ValueError('{} is not equivalent to {}'.format(unit, given_unit))
        self._uncertainty = u.Quantity(value, unit).to(given_unit)

    @property
    def reference(self):
        return self._reference

    @reference.setter
    def reference(self, value):
        if value is None:
            value = ""
        if not isinstance(value, str):
            raise TypeError('reference must be a string')
        self._reference = value

    @property
    def notes(self):
        return self._notes

    @notes.setter
    def notes(self, value):
        if value is None:
            value = ""
        if not isinstance(value, str):
            raise TypeError('notes must be a string')
        self._notes = value

    def __repr__(self):
        if not np.isnan(self.value):
            return super().__repr__()
        else:
            return ""

    def __str__(self):
        if not np.isnan(self.value):
            out = "{}:\n    {:.5g} +/- {:.5g} {}\n    Reference: {}, {}\n".format(
                self.name, self.value, self.uncertainty.value, self.unit.__str__(), self.reference, self.notes)
            return out
        else:
            return ""


class BaseBody():

    @property
    def albedo(self):
        return self._albedo

    @albedo.setter
    def albedo(self, value):
        if isinstance(value, PhysicalData):
            self._albedo = value
        else:
            self._albedo = PhysicalData('Albedo', value)

    @property
    def H(self):
        return self._H

    @H.setter
    def H(self, value):
        if isinstance(value, PhysicalData):
            self._H = value
        else:
            self._H = PhysicalData('Absolute Magnitude', value, unit=u.mag)
        if not np.isnan(self.H):  # remove this line for v1.0
            self._shared_with['ephem']['H'] = self.H  # remove this line for v1.0

    @property
    def G(self):
        return self._G

    @G.setter
    def G(self, value):
        if isinstance(value, PhysicalData):
            self._G = value
        else:
            self._G = PhysicalData('Phase Slope', value)
        if not np.isnan(self._G):  # remove this line for v1.0
            self._shared_with['ephem']['G'] = self.G  # remove this line for v1.0

    @property
    def diameter(self):
        return self._diameter

    @diameter.setter
    def diameter(self, value):
        if isinstance(value, PhysicalData):
            self._diameter = value
        else:
            self._diameter = PhysicalData('Diameter', value, unit=u.km)
        self._shared_with['ephem']['radius'] = self.radius

    @property
    def radius(self):
        return PhysicalData('Radius', self.diameter/2.0, self.diameter.uncertainty/2,
                            self.diameter.reference, self.diameter.notes, unit=u.km)

    @radius.setter
    def radius(self, value):
        if isinstance(value, PhysicalData):
            self.diameter = value*2.0
            self.diameter.name = 'Diameter'
        else:
            self.diameter = PhysicalData('Diameter', value*2.0, unit=u.km)

    @property
    def density(self):
        return self._density

    @density.setter
    def density(self, value):
        if isinstance(value, PhysicalData):
            self._density = value
        else:
            self._density = PhysicalData('Density', value, unit=u.g/u.cm**3)

    @property
    def GM(self):
        return self._GM

    @GM.setter
    def GM(self, value):
        if isinstance(value, PhysicalData):
            self._GM = value
        else:
            self._GM = PhysicalData('Standard Gravitational Parameter', value, unit=u.km**3/u.s**2)

    @property
    def mass(self):
        return PhysicalData('Mass', self.GM/const.G, self.GM.uncertainty/const.G,
                            self.GM.reference, self.GM.notes, unit=u.kg)

    @property
    def rotation(self):
        return self._rotation

    @rotation.setter
    def rotation(self, value):
        if isinstance(value, PhysicalData):
            self._rotation = value
        else:
            self._rotation = PhysicalData('Rotation', value, unit=u.h)

    @property
    def pole(self):
        return self._pole

    @pole.setter
    def pole(self, value):
        if value is None:
            self._pole = SkyCoord(np.nan, np.nan, unit=(u.deg, u.deg))
        else:
            self._pole = SkyCoord(value, unit=(u.hourangle, u.deg))
        self._pole.ra.uncertainty = Longitude(0*u.hourangle)
        self._pole.dec.uncertainty = Latitude(0*u.deg)
        self._pole.reference = "User"
        self._pole.notes = ""

    @property
    def BV(self):
        return self._BV

    @BV.setter
    def BV(self, value):
        if isinstance(value, PhysicalData):
            self._BV = value
        else:
            self._BV = PhysicalData('B-V color', value)

    @property
    def UB(self):
        return self._UB

    @UB.setter
    def UB(self, value):
        if isinstance(value, PhysicalData):
            self._UB = value
        else:
            self._UB = PhysicalData('U-B color', value)

    @property
    def spkid(self):
        return self._shared_with['ephem'].get('spkid')

    @spkid.setter
    def spkid(self, value):
        spkid = str(int(value))
        self._shared_with['ephem']['spkid'] = spkid

    @property
    def ephem(self):
        try:
            return self._ephem
        except AttributeError:
            raise AttributeError('{} object does not have ephemeris.'.format(self.__class__.__name__))

    @ephem.setter
    def ephem(self, value):
        allowed_types = [EphemPlanete, EphemKernel, EphemJPL, EphemHorizons]
        if type(value) not in allowed_types:
            if isinstance(value, str) and value.lower() == 'horizons':
                value = EphemHorizons(name=self._search_name)
            elif isinstance(value, (list, str)):
                value = EphemKernel(kernels=value, spkid=self.spkid)
            else:
                raise ValueError('Cannot set "ephem" with {}. Allowed types are: {}'.format(type(value), allowed_types))
        if 'spkid' not in self._shared_with['ephem'] or self._shared_with['ephem']['spkid'] is None:
            if hasattr(value, '_spkid'):
                self._shared_with['ephem']['spkid'] = value.spkid
            else:
                raise AttributeError('spkid is not defined in {} or {}'.format(self.__class__.__name__, value.__class__.__name__))
        if hasattr(self, '_ephem'):
            self._ephem._shared_with['body'] = {}
        self._ephem = value
        spkval = value.spkid
        self._ephem._shared_with['body'] = self._shared_with['ephem']
        spknewval = value.spkid
        if spkval != spknewval:
            warnings.warn('spkid is different in {0} ({1}) and {2} ({3}). {0}\'s spkid will have higher priority'.format(
                self.__class__.__name__, spknewval, value.__class__.__name__, spkval))

    @property
    def _search_name(self):
        if self.orbit_class in ['Natural Satellite', 'Planet']:
            return str(self.spkid or self.shortname)
        elif hasattr(self, '_des_name'):
            return self._des_name
        else:
            return self.name

    @property
    def _id_type(self):
        if self.orbit_class in ['Natural Satellite', 'Planet']:
            return 'majorbody'
        else:
            return 'smallbody'
