import astropy.units as u
import numpy as np


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
        cls = super().__new__(cls, value=value, unit=unit)
        cls.name = name
        cls.uncertainty = uncertainty
        cls.reference = reference
        cls.notes = notes
        return cls

    @property
    def uncertainty(self):
        return self._uncertainty

    @uncertainty.setter
    def uncertainty(self, value):
        unit = self.unit
        if isinstance(value, u.core.CompositeUnit):
            unit = u.dimensionless_unscaled
            for i in np.arange(len(value.bases)):
                unit = unit*(value.bases[i]**value.powers[i])
            value = value.scale
        elif value is None:
            value = 0.0
        if not unit.is_equivalent(self.unit):
            raise u.UnitError('Uncertainty unit must be equivalent to {}'.format(self.unit.name))
        self._uncertainty = u.Quantity(value, unit)

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
