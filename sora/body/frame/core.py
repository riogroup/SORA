import astropy.units as u
import numpy as np
from astropy.coordinates import BaseCoordinateFrame, frame_transform_graph, AffineTransform, ICRS, SkyCoord, \
    CartesianRepresentation, Angle, SphericalRepresentation
from astropy.coordinates.attributes import CoordinateAttribute, TimeAttribute, QuantityAttribute, Attribute
from astropy.time import Time
from .meta import Precession

__all__ = ['PlanetocentricFrame']


class PlanetocentricFrame(BaseCoordinateFrame):
    """
    https://docs.astropy.org/en/stable/coordinates/frames.html

    Parameters
    ----------
    epoch : `str`, `astropy.time.Time`
        Reference epoch of the given parameters, in TDB.

    pole : `str`, `astropy.coordinates.SkyCoord`
        ICRS coordinates of the pole at reference epoch.
        If string, it must be 'hh.hhhh +dd.ddd' or 'hh hh hh.hhh +dd dd dd.ddd',
        in hourangle and deg.

    alphap : `float`, `astropy.units.Quantity`
        Rate at which the right ascension of the pole changes, in deg/century

    extra_alpha: `Precession`
        A Precession object that accounts for the precession of the pole in right ascension

    deltap : `float`, `astropy.units.Quantity`
        Rate at which the declination of the pole changes, in deg/century

    extra_delta : `Precession`
        A Precession object that accounts for the precession of the pole in right ascension

    prime_angle : `float`, `astropy.units.Quantity`
        The angle of the prime meridian at reference epoch, in deg

    rotation_velocity : `float`, `astropy.units.Quantity`
        The rotation velocity of the body, in deg/day

    extra_w : `Precession`
        A Precession object that accounts for the precession of the pole in right ascension

    right_hand : `bool`
        States if the orientation of the longitude of the body is counterclockwise.
        In the Solar System, Earth, the Moon and the Sun have longitude counterclockwise.
        The asteroids, planets and satellites must define to have an increasing longitude
        from the Earth's point-of-view.

    reference : `str`
        Includes a reference or citation for the given parameters

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
        """

        Parameters
        ----------
        epoch : `str`, `astropy.time.Time`
            Time to which rotate the frame.

        Returns
        -------
        pole : `astropy.coordinates.SkyCoord`
            The pole of the object at given epoch.

        W : `astropy.coordinates.Angle`
            The location of the prime meridian relative to the ascending node of the body's equator
            at given epoch.

        """
        dt = Time(epoch) - self.epoch
        W = self.prime_angle + self.rotation_velocity * dt + self.extra_w.compute_at(dt)
        W = Angle(W).wrap_at(360 * u.deg)
        new_ra = self.pole.ra + self.alphap * dt + self.extra_alpha.compute_at(dt)
        new_dec = self.pole.dec + self.deltap * dt + self.extra_delta.compute_at(dt)
        pole = SkyCoord(new_ra, new_dec, frame='icrs')
        return pole, W

    def frame_at(self, epoch):
        """

        Parameters
        ----------
        epoch : `str`, `astropy.time.Time`
            Time to which rotate the frame.

        Returns
        -------
        frame : `PlanetocentricFrame`
            A new PlanetocentricFrame with the parameters at given epoch.

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
                  "    Epoch: {}".format(self.epoch.__str__()),
                  "    alpha_pole = {} {:+f}*T {}".format(self.pole.ra.value, self.alphap.value*100,
                                                          ''.join(self.extra_alpha.__str__().split('\n'))),
                  "    delta_pole = {} {:+f}*T {}".format(self.pole.dec.value, self.deltap.value*100,
                                                          ''.join(self.extra_delta.__str__().split('\n'))),
                  "    W = {} {:+f}*d {}".format(self.prime_angle.value, self.rotation_velocity.value,
                                                 ''.join(self.extra_w.__str__().split('\n'))),
                  "    Reference: {}".format(self.reference)]
        return '\n'.join(string)


def get_matrix_vectors(planetocentric_frame, inverse=False):
    """

    Parameters
    ----------
    planetocentric_frame : `PlanetocentricFrame`
        The PlanetocentricFrame object to convert from the ICRS

    inverse : `bool`
        If the parameters are to be calculated to the ICRS

    Returns
    -------
    : `numpy.array`
        A matrix to convert between orientations

    : `CartesianRepresentation`
        A vector to convert between origins

    """
    from astropy.coordinates.matrix_utilities import rotation_matrix, matrix_product, matrix_transpose
    offset = CartesianRepresentation(0 * u.km, 0 * u.km, 0 * u.km)
    pole, W = planetocentric_frame.orientation_at(planetocentric_frame.epoch)
    if planetocentric_frame.right_hand:
        m1 = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    else:
        m1 = np.array([[1, 0, 0], [0, -1, 0], [0, 0, 1]])
    rz2 = rotation_matrix(W, axis='z')
    rx = rotation_matrix(90 * u.deg - pole.dec, axis='x')
    rz1 = rotation_matrix(pole.ra + 90 * u.deg, axis='z')

    A = matrix_product(m1, rz2, rx, rz1)
    if inverse:
        A = matrix_transpose(A)
    return A, offset


@frame_transform_graph.transform(AffineTransform, ICRS, PlanetocentricFrame)
def icrs_to_planetocentric(icrs_coord, planetocentric_frame):
    return get_matrix_vectors(planetocentric_frame)


@frame_transform_graph.transform(AffineTransform, PlanetocentricFrame, ICRS)
def planetocentric_to_icrs(planetocentric_coord, icrs_frame):
    return get_matrix_vectors(planetocentric_coord, inverse=True)
