# sora/rings/core.py
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
        # --- attach body if provided ---
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

        # If no explicit pole is given, inherit from body if possible
        if pole is None and self.body is not None:
            body_pole = getattr(self.body, "pole", None)
            if body_pole is not None and not np.isnan(body_pole.ra.deg):
                pole = body_pole

        # Build RingGeometry
        if pole is None:
            self.geometry = RingGeometry(pole_ra=None, pole_dec=None)
        else:
            pole = SkyCoord(pole)
            self.geometry = RingGeometry(pole_ra=pole.ra.deg, pole_dec=pole.dec.deg)



        # --- physical properties (BaseRing init) ---
        super().__init__(
            radius=kwargs.get('radius'),
            radial_width=kwargs.get('radial_width'),
            normal_opacity=kwargs.get('normal_opacity'),
            normal_optical_depth=kwargs.get('normal_optical_depth'),
            equivalent_depth=kwargs.get('equivalent_depth'),
            equivalent_width=kwargs.get('equivalent_width'),
            eccentricity=kwargs.get('eccentricity'),
        )

        # --- attach uncertainties ---
        self._radius.uncertainty = kwargs.get('radius_err', 0.0)
        self._radial_width.uncertainty = kwargs.get('radial_width_err', 0.0)
        self._normal_opacity.uncertainty = kwargs.get('normal_opacity_err', 0.0)
        self._normal_optical_depth.uncertainty = kwargs.get('normal_optical_depth_err', 0.0)
        self._equivalent_depth.uncertainty = kwargs.get('equivalent_depth_err', 0.0)
        self._equivalent_width.uncertainty = kwargs.get('equivalent_width_err', 0.0)
        self._eccentricity.uncertainty = kwargs.get('eccentricity_err', 0.0)

    # ------------------------------------------------------------------
    def get_ring_orientation(self, time, observer="geocenter"):
        return self.geometry.orientation(self.ephem, time, observer)

    # ------------------------------------------------------------------
    def to_ring_plane(self, f, g, time, center_f=0, center_g=0):
        """
        Convert sky-plane coordinates (f, g) to ring-plane (x, y),
        computing internally the projection coefficients.
        """
        pos = self.ephem.get_position(time)
        P, B = self.get_ring_orientation(time)
        earth_pole = SkyCoord('12h00m00s +90d00m00s')

        coef, coef_polo = calc_coef_projecao(pos, self.geometry.pole, B, P, earth_pole)
        x, y = project_to_ring_plane(f, g, coef, coef_polo, ksi_0=center_f, eta_0=center_g)
        return x, y

    # ------------------------------------------------------------------
    def __str__(self):
        out = []
        out.append(f"Ring ID: {self.ring_id}\n")
        out.append(str(self.geometry))
        out.append(self._radius.__str__() + "\n")
        out.append(self._normal_opacity.__str__() + "\n")
        out.append(self._normal_optical_depth.__str__() + "\n")
        out.append(self._radial_width.__str__() + "\n")
        out.append(self._eccentricity.__str__() + "\n")
        out.append(self._equivalent_width.__str__() + "\n")
        out.append(self._equivalent_depth.__str__() + "\n")
        return ''.join(out)
