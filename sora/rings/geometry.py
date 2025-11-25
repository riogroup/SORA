# sora/rings/geometry.py
import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time

__all__ = ["RingGeometry"]


class RingGeometry:
    """
    RingGeometry — geometric and orientation parameters of a ring.

    Handles pole orientation, projected opening angle (B), and position
    angle of the ring plane (P).

    Parameters
    ----------
    pole_ra : float or None
        Right Ascension of the ring pole (degrees or hourangle).
    pole_dec : float or None
        Declination of the ring pole (degrees).
    epoch : `astropy.time.Time`, optional
        Epoch for which the orientation is defined.
    reference : str, optional
        Reference or publication of the adopted pole.

    Notes
    -----
    - This class does not include the physical properties of the ring
      (radius, opacity, etc.), which are handled by `BaseRing`.
    - It is designed to be embedded in `Ring` objects.
    """

    def __init__(self, pole_ra=None, pole_dec=None, epoch=None, reference="User"):
        self.epoch = epoch
        if pole_ra is None or pole_dec is None:
            self.pole = SkyCoord(np.nan, np.nan, unit=(u.deg, u.deg))
        else:
            self.pole = SkyCoord(ra=pole_ra*u.deg, dec=pole_dec*u.deg, frame="icrs")


    def orientation(self, ephem, time, observer="geocenter"):
        time = Time(time)
        pos = ephem.get_position(time, observer=observer)

        pole = self.pole    

        P = pos.position_angle(pole).to(u.deg)
        B = np.arcsin(
            -(np.sin(pole.dec)*np.sin(pos.dec) +
            np.cos(pole.dec)*np.cos(pos.dec)*np.cos(pole.ra - pos.ra))
        ).to(u.deg)

        return P, B

    def __str__(self):
        out = "Ring Geometry:\n"
        out += f"  Pole (RA, Dec): {self.pole.ra.deg:.3f}°, {self.pole.dec.deg:.3f}°\n"
        out += f"  Reference: {self.reference}\n"
        return out
