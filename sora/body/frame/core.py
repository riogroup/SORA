import astropy.units as u
import numpy as np
from astropy.coordinates import BaseCoordinateFrame, frame_transform_graph, AffineTransform, ICRS, SkyCoord, \
    CartesianRepresentation, Angle, SphericalRepresentation
from astropy.coordinates.attributes import CoordinateAttribute, TimeAttribute, QuantityAttribute, Attribute
from astropy.time import Time
from .meta import Precession

__all__ = ['PlanetocentricFrame']


class PlanetocentricFrame(BaseCoordinateFrame):
    """Represent a planetocentric coordinate frame.

    This frame stores the pole orientation, prime meridian angle, rotation
    velocity, and optional precession terms needed to transform between ICRS and
    body-fixed coordinates.

    Parameters
    ----------
    epoch : `str`, `astropy.time.Time`
        Reference epoch of the given parameters.

    pole : `str`, `astropy.coordinates.SkyCoord`
        ICRS coordinates of the pole at reference epoch.
        If a string is given, it must be in the format
        ``'hh.hhhh +dd.ddd'`` or ``'hh hh hh.hhh +dd dd dd.ddd'``, in
        hourangle and degrees.

    alphap : `float`, `astropy.units.Quantity`
        Rate at which the right ascension of the pole changes, in degrees per
        century.

    extra_alpha : `Precession`
        Precession terms for the pole right ascension.

    deltap : `float`, `astropy.units.Quantity`
        Rate at which the declination of the pole changes, in degrees per
        century.

    extra_delta : `Precession`
        Precession terms for the pole declination.

    prime_angle : `float`, `astropy.units.Quantity`
        Angle of the prime meridian at reference epoch, in degrees.

    rotation_velocity : `float`, `astropy.units.Quantity`
        Rotation velocity of the body, in degrees per day.

    extra_w : `Precession`
        Precession terms for the prime meridian angle.

    right_hand : `bool`
        If True, the longitude orientation follows the right-hand rule. In the
        Solar System, Earth, the Moon, and the Sun use counterclockwise
        longitudes. Asteroids, planets, and satellites must be defined to have
        increasing longitude from the Earth's point of view.

    reference : `str`
        Reference or citation for the given parameters.

    """
    # Specify how coordinate values are represented when outputted
    default_representation = SphericalRepresentation

    epoch = TimeAttribute(default='J2000')
    pole = CoordinateAttribute(default=SkyCoord("00 90", unit=('hourangle', 'deg')), frame=ICRS)
    alphap = QuantityAttribute(default=0 * u.deg / u.year)
    extra_alpha = Attribute(default=Precession())
    deltap = QuantityAttribute(default=0 * u.deg / u.year)
    extra_delta = Attribute(default=Precession())
    prime_angle = QuantityAttribute(default=0 * u.deg)
    rotation_velocity = QuantityAttribute(default=0 * u.deg / u.day)
    extra_w = Attribute(default=Precession())
    right_hand = Attribute(default=False)
    reference = Attribute(default="User")

    def __init__(self, *args, **kwargs):
        pole = kwargs.get('pole')
        if isinstance(pole, str):
            kwargs['pole'] = SkyCoord(pole, unit=(u.hourangle, u.deg), frame='icrs')
        kwargs['alphap'] = u.Quantity(kwargs.get('alphap', 0), unit=u.deg / u.year) / 100
        kwargs['deltap'] = u.Quantity(kwargs.get('deltap', 0), unit=u.deg / u.year) / 100
        kwargs['prime_angle'] = Angle(u.Quantity(kwargs.get('prime_angle', 0), unit=u.deg)).wrap_at(360 * u.deg)
        kwargs['rotation_velocity'] = u.Quantity(kwargs.get('rotation_velocity', 0), unit=u.deg / u.day)
        kwargs['extra_alpha'] = Precession(kwargs.get('extra_alpha', 0), func='sin', multiplier='T')
        kwargs['extra_delta'] = Precession(kwargs.get('extra_delta', 0), func='cos', multiplier='T')
        kwargs['extra_w'] = Precession(kwargs.get('extra_w', 0), func='sin', multiplier='T')
        super().__init__(*args, **kwargs)

    def orientation_at(self, epoch):
        """Return the pole and prime meridian orientation at an epoch.

        Parameters
        ----------
        epoch : `str`, `astropy.time.Time`
            Time at which to evaluate the frame orientation. If a string is
            given, the scale is UTC.

        Returns
        -------
        pole : `astropy.coordinates.SkyCoord`
            Pole of the object at the given epoch.

        W : `astropy.coordinates.Angle`
            Location of the prime meridian relative to the ascending node of the
            body's equator at the given epoch.

        """
        dt = Time(epoch) - self.epoch
        W = self.prime_angle + self.rotation_velocity * dt + self.extra_w.compute_at(dt)
        W = Angle(W).wrap_at(360 * u.deg)
        new_ra = self.pole.ra + self.alphap * dt + self.extra_alpha.compute_at(dt)
        new_dec = self.pole.dec + self.deltap * dt + self.extra_delta.compute_at(dt)
        pole = SkyCoord(new_ra, new_dec, frame='icrs')
        return pole, W

    def frame_at(self, epoch):
        """Return a frame propagated to a new epoch.

        Parameters
        ----------
        epoch : `str`, `astropy.time.Time`
            Time at which to evaluate the frame. If a string is given, the scale
            is UTC.

        Returns
        -------
        frame : `PlanetocentricFrame`
            New frame with parameters propagated to the given epoch.

        """
        dt = Time(epoch) - self.epoch
        pole, W = self.orientation_at(epoch=epoch)
        alpha = pole.ra - self.extra_alpha.compute_at(dt)
        delta = pole.dec - self.extra_delta.compute_at(dt)
        new_pole = SkyCoord(alpha, delta)
        W = Angle(W - self.extra_w.compute_at(dt)).wrap_at(360 * u.deg)
        extra_alpha = self.extra_alpha.params_at(dt)
        extra_delta = self.extra_delta.params_at(dt)
        extra_w = self.extra_w.params_at(dt)
        new_frame = PlanetocentricFrame(epoch=epoch, pole=new_pole, alphap=self.alphap, deltap=self.deltap, prime_angle=W,
                                        rotation_velocity=self.rotation_velocity, right_hand=self.right_hand,
                                        reference=self.reference, extra_w=extra_w, extra_alpha=extra_alpha,
                                        extra_delta=extra_delta)
        return new_frame

    def __str__(self):
        string = ["PlanetocentricFrame:",
                  "    Epoch: {} {}".format(self.epoch.__str__(), self.epoch.scale),
                  "    alpha_pole = {} {:+f}*T {}".format(self.pole.ra.value, self.alphap.value*100,
                                                          ''.join(self.extra_alpha.__str__().split('\n'))),
                  "    delta_pole = {} {:+f}*T {}".format(self.pole.dec.value, self.deltap.value*100,
                                                          ''.join(self.extra_delta.__str__().split('\n'))),
                  "    W = {} {:+f}*d {}".format(self.prime_angle.value, self.rotation_velocity.value,
                                                 ''.join(self.extra_w.__str__().split('\n'))),
                  "    Reference: {}".format(self.reference)]
        return '\n'.join(string)


def get_matrix_vectors(planetocentric_frame, inverse=False):
    """Return the transformation matrix and offset vector for a frame.

    Parameters
    ----------
    planetocentric_frame : `PlanetocentricFrame`
        Planetocentric frame used to calculate the transform from ICRS.

    inverse : `bool`
        If True, return the inverse transform to ICRS.

    Returns
    -------
    matrix : `numpy.ndarray`
        Matrix used to convert between orientations.

    offset : `astropy.coordinates.CartesianRepresentation`
        Vector used to convert between origins.

    """
    from astropy.coordinates.matrix_utilities import rotation_matrix, matrix_transpose
    offset = CartesianRepresentation(0 * u.km, 0 * u.km, 0 * u.km)
    pole, W = planetocentric_frame.orientation_at(planetocentric_frame.epoch)
    if planetocentric_frame.right_hand:
        m1 = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    else:
        m1 = np.array([[1, 0, 0], [0, -1, 0], [0, 0, 1]])
    rz2 = rotation_matrix(W, axis='z')
    rx = rotation_matrix(90 * u.deg - pole.dec, axis='x')
    rz1 = rotation_matrix(pole.ra + 90 * u.deg, axis='z')

    A = m1 @ rz2 @ rx @ rz1
    if inverse:
        A = matrix_transpose(A)
    return A, offset


@frame_transform_graph.transform(AffineTransform, ICRS, PlanetocentricFrame)
def icrs_to_planetocentric(icrs_coord, planetocentric_frame):
    """Return the affine transform from ICRS to a planetocentric frame."""
    return get_matrix_vectors(planetocentric_frame)


@frame_transform_graph.transform(AffineTransform, PlanetocentricFrame, ICRS)
def planetocentric_to_icrs(planetocentric_coord, icrs_frame):
    """Return the affine transform from a planetocentric frame to ICRS."""
    return get_matrix_vectors(planetocentric_coord, inverse=True)
