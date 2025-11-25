# sora/rings/core.py


# TODO (REFORMAT-PHYSICALDATA):
# The classes storing ring physical properties (e.g., PhysicalData-derived fields)
# generate excessive blank lines when optional attributes (error, reference, notes)
# are None. Their __str__() methods should be refactored to:
#
# 1. Print only fields that have meaningful values.
# 2. Avoid empty lines caused by missing metadata.
# 3. Remove trailing commas or partially filled lines (e.g., "Reference: User, ").
# 4. Return an empty string when the property value itself is None.
#
# This affects the formatted output of Ring.__str__() because each PhysicalData
# block inserts extra newlines even when nearly all fields are undefined.
#
# Proposed fix:
# - Update PhysicalData.__str__() to conditionally assemble lines:
#     - Always show the property name and value.
#     - Only show error, reference, notes when they are not None/empty.
#     - Strip trailing whitespace and collapse empty lines.
#
# This will eliminate large blank gaps in print(Ring) and produce a clean,
# compact summary of the ring's parameters.


import warnings
import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time

from sora.config import input_tests
from sora.body.meta import PhysicalData
from sora.body import Body
from sora.ephem.meta import BaseEphem

from .meta import BaseRing
from .geometry import RingGeometry
from .utils import calc_coef_projecao, project_to_ring_plane


__all__ = ["Ring"]


class Ring(BaseRing):
    """
    Ring — container for the physical and geometric information of a ring.

    Compatible with the previous implementation, but internally structured
    around BaseRing (physical parameters) and RingGeometry (orientation).
    """

    def __init__(self, **kwargs):
        self.body = kwargs.pop("body", None)
        if self.body is not None:
            self.ephem = self.body.ephem
        else:
            self.ephem = None

        allowed_kwargs = [
            'ring_id',
            'radius', 'radius_err',
            'eccentricity', 'eccentricity_err',
            'pole_orientation',
            'normal_opacity', 'normal_opacity_err',
            'normal_optical_depth', 'normal_optical_depth_err',
            'radial_width', 'radial_width_err',
            'equivalent_depth', 'equivalent_depth_err',
            'equivalent_width', 'equivalent_width_err'
        ]

        kwargs.pop("body", None)
        input_tests.check_kwargs(kwargs, allowed_kwargs=allowed_kwargs)

        self.ring_id = kwargs.get('ring_id', 'Unknown')

        pole = kwargs.get("pole_orientation", None)

        if pole is None and self.body is not None:
            body_pole = getattr(self.body, "pole", None)
            if body_pole is not None and not np.isnan(body_pole.ra.deg):
                pole = body_pole

        if pole is None:
            self.geometry = RingGeometry(pole_ra=None, pole_dec=None)
        else:
            pole = SkyCoord(pole)
            self.geometry = RingGeometry(pole_ra=pole.ra.deg, pole_dec=pole.dec.deg)


        super().__init__(
            radius=kwargs.get('radius'),
            radial_width=kwargs.get('radial_width'),
            normal_opacity=kwargs.get('normal_opacity'),
            normal_optical_depth=kwargs.get('normal_optical_depth'),
            equivalent_depth=kwargs.get('equivalent_depth'),
            equivalent_width=kwargs.get('equivalent_width'),
            eccentricity=kwargs.get('eccentricity'),
        )

        self._radius.uncertainty = kwargs.get('radius_err', 0.0)
        self._radial_width.uncertainty = kwargs.get('radial_width_err', 0.0)
        self._normal_opacity.uncertainty = kwargs.get('normal_opacity_err', 0.0)
        self._normal_optical_depth.uncertainty = kwargs.get('normal_optical_depth_err', 0.0)
        self._equivalent_depth.uncertainty = kwargs.get('equivalent_depth_err', 0.0)
        self._equivalent_width.uncertainty = kwargs.get('equivalent_width_err', 0.0)
        self._eccentricity.uncertainty = kwargs.get('eccentricity_err', 0.0)

    def get_ring_orientation(self, time, observer="geocenter"):
        """
        Return the instantaneous ring orientation as seen by the specified observer.

        This method computes the projected position angle (P) and opening angle (B)
        of the ring plane on the sky at a given epoch. The computation is delegated
        to the `RingGeometry.orientation()` method.

        Parameters
        ----------
        time : str or astropy.time.Time
            Epoch at which the ring orientation is evaluated.
        observer : str or sora.Observer, optional
            Observer location. Defaults to 'geocenter'.

        Returns
        -------
        P, B : astropy.units.Quantity
            - P : position angle of the ring plane (degrees), measured east of north.
            - B : opening angle of the ring (degrees); B = 0° corresponds to edge-on.

        """
        return self.geometry.orientation(self.ephem, time, observer)

    def to_ring_plane(self, f, g, time, center_f=0, center_g=0, observer="geocenter"):
        """
        Convert sky-plane coordinates (f, g) into the ring-plane coordinates (x, y).

        This method transforms observer-centered sky-plane positions into the 
        equatorial plane of the ring using the instantaneous ring orientation. 
        Internally, it:

        1. Obtains the apparent position of the body at the given epoch (time).
        2. Computes the ring opening angle (B) and position angle (P).
        3. Applies the projection coefficients to map (f, g) → (x, y).

        Parameters
        ----------
        f, g : float or array_like
            Sky-plane coordinates in km (positive f = celestial east, positive g = celestial north).
        time : str or astropy.time.Time
            Epoch at which the projection is computed.
        center_f, center_g : float, optional
            Offsets of the body's center in the (f, g) plane, in km. 
            Defaults to 0 for both.
        observer : str or sora.Observer, optional
            Observer location. Defaults to 'geocenter'.

        Returns
        -------
        x, y : ndarray
            Coordinates in the ring-plane (equatorial plane of the ring), in km.

        Notes
        -----
        - The result depends on the instantaneous ring orientation relative to the
        observer, meaning (x, y) are time-dependent quantities.
        """
        pos = self.ephem.get_position(time, observer=observer)
        P, B = self.geometry.orientation(self.ephem, time, observer)
        return self.geometry.to_ring_plane(pos, f, g, P, B, center_f, center_g)
    
    def __str__(self):
        out = [f"Ring ID: {self.ring_id}\n"]
        out.append(str(self.geometry))

        props = [
            self._radius,
            self._normal_opacity,
            self._normal_optical_depth,
            self._radial_width,
            self._eccentricity,
            self._equivalent_width,
            self._equivalent_depth,
        ]

        for prop in props:
            if prop is not None and prop.value is not None:
                out.append(str(prop) + "\n")

        return "".join(out)

