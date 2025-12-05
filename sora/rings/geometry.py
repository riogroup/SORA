# sora/rings/geometry.py
import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord, Longitude, Latitude
from astropy.time import Time
from .utils import calc_coef_projecao, project_to_ring_plane

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

    def __init__(self, pole_ra=None, pole_dec=None, epoch=None):
        """
        Initialize the geometric description of a ring.

        Parameters
        ----------
        pole_ra : float or None
            Right ascension of the ring pole, in degrees. If None, the pole
            is initialized as undefined (NaN).
        pole_dec : float or None
            Declination of the ring pole, in degrees. If None, the pole
            is initialized as undefined (NaN).
        epoch : astropy.time.Time or None, optional
            Epoch associated with the pole orientation (if relevant for
            external references). Currently not used in the projection
            routines, but preserved for metadata consistency.

        Notes
        -----
        - The ring pole is stored internally as an ICRS `SkyCoord`.
        - A pole defined as NaN indicates that the ring orientation is
        not yet specified and must be provided before geometric
        projections or orientation calculations can be performed.
        """
        self.epoch = epoch
        if pole_ra is None or pole_dec is None:
            self._pole = SkyCoord(np.nan, np.nan, unit=(u.deg, u.deg))
        else:
            self._pole = SkyCoord(ra=pole_ra*u.deg, dec=pole_dec*u.deg)

    @property
    def pole_orientation(self):
        """Return the ring pole as a SkyCoord."""
        return self._pole

    @pole_orientation.setter
    def pole_orientation(self, value):
        """Set the ring pole (SkyCoord or parseable string)."""
        if value is None:
            self._pole = SkyCoord(np.nan, np.nan, unit=(u.deg, u.deg))
        else:
            self._pole = SkyCoord(value, unit=(u.hourangle, u.deg))
        self._pole.reference = "User"
    
    def orientation(self, ephem, time, observer="geocenter"):
        """
        Compute the apparent orientation of the ring as seen by the observer.

        This method returns the two standard geometric angles that define the
        projected appearance of a ring on the sky:

        - P : position angle of the ring plane, measured east of celestial north.
        - B : opening angle of the ring; B = 0° corresponds to an edge-on view.

        Both quantities depend on the apparent position of the central body
        (from `ephem.get_position`) and on the fixed pole of the ring.

        Parameters
        ----------
        ephem : Ephem instance
            Ephemeris object associated with the central body.
        time : str or astropy.time.Time
            Epoch at which the ring orientation is evaluated.
        observer : str or sora.Observer, optional
            Observer code or object. Default is 'geocenter'.

        Returns
        -------
        P, B : astropy.units.Quantity
            P : position angle of the ring plane in degrees.
            B : opening angle of the ring plane in degrees.
        """

        time = Time(time)
        pos = ephem.get_position(time, observer=observer)

        pole_orientation = self.pole_orientation    

        P = pos.position_angle(pole_orientation).to(u.deg)
        B = np.arcsin(
            -(np.sin(pole_orientation.dec)*np.sin(pos.dec) +
            np.cos(pole_orientation.dec)*np.cos(pos.dec)*np.cos(pole_orientation.ra - pos.ra))
        ).to(u.deg)

        return P, B
    
    def to_ring_plane(self, pos, f, g, P, B, center_f=0, center_g=0, earth_pole = SkyCoord('12h00m00s +90d00m00s')):
        """
        Convert sky-plane coordinates (f, g) to ring-plane (x, y).

        Parameters
        ----------
        pos : SkyCoord
            Apparent position of the body (output from ephem.get_position()).
        f, g : array_like
            Sky-plane coords (km)
        P, B : Quantity
            Position angle, opening angle (degrees)
        center_f, center_g : float
            Offsets of body center in km.
        """

        coef, coef_polo = calc_coef_projecao(
            ephem=pos,              
            pole_coord=self.pole_orientation,
            B=B,
            P=P,
            earth_pole=earth_pole
        )

        x, y = project_to_ring_plane(
            f, g,
            coef, coef_polo,
            ksi_0=center_f,
            eta_0=center_g
        )

        return x, y

    def __str__(self):
        if not np.isnan(self._pole.ra):
            out = ('Pole orientation\n    RA:{}\n    DEC:{}\n    Reference: {}\n'.format(
                self._pole.ra.__str__(), 
                self._pole.dec.__str__(),
                self._pole.reference.__str__()))
            return ''.join(out)
        return ''
        
