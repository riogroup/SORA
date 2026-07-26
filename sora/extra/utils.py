import numpy as np

__all__ = ['get_ellipse_points']


def get_ellipse_points(theta, equatorial_radius, oblateness=0.0, center_f=0.0, center_g=0.0,
                       position_angle=0.0):
    """Calculate points on an ellipse for the given input parameters.

    Parameters
    ----------
    theta : `float`, array-like
        Angular coordinate, in radians, where the ellipse points are
        calculated.

    equatorial_radius : `float`, `int`
        Semi-major axis of the ellipse, in km.

    oblateness : `float`, `int`, optional
        Oblateness of the ellipse.

    center_f : `float`, `int`, optional
        Coordinate of the ellipse center in the f direction, in km.

    center_g : `float`, `int`, optional
        Coordinate of the ellipse center in the g direction, in km.

    position_angle : `float`, `int`, optional
        The pole position angle of the ellipse in degrees.
        Zero is in the North direction ('g-positive'). Positive clockwise.

    Returns
    -------
    x_model : `float`, `numpy.array`
        Cartesian x-component, in km.

    y_model : `float`, `numpy.array`
        Cartesian y-component, in km.

    r_model : `float`, `numpy.array`
        Radial distance, in km.

    theta : `float`, array-like
        Input angular coordinate, in radians.
    """
    a = equatorial_radius
    b = equatorial_radius - equatorial_radius * oblateness
    phi = position_angle * (np.pi / 180.0)
    ang = theta + phi
    r_model = (a * b) / np.sqrt((a * np.sin(ang)) ** 2 + (b * np.cos(ang)) ** 2)
    x_model = r_model * np.cos(theta) + center_f
    y_model = r_model * np.sin(theta) + center_g
    return x_model, y_model, r_model, theta
