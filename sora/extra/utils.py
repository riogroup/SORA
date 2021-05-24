import numpy as np

__all__ = ['get_ellipse_points']


def get_ellipse_points(theta, equatorial_radius, oblateness=0.0, center_f=0.0, center_g=0.0,
                       position_angle=0.0):
    """Get points for the ellipse with the given input parameters.

    Parameters
    ----------
    theta : `float` array
        Angular coordinate, in degrees, to return the ellipse points.

    equatorial_radius : `float`, `int`
        Semi-major axis of the ellipse.

    oblateness : `float`, `int`, default=0
        Oblateness of the ellipse.

    center_f : `float`, `int`, default=0
        Coordinate of the ellipse (abscissa).

    center_g : `float`, `int`, default=0
        Coordinate of the ellipse (ordinate).

    position_angle : `float`, `int`, default=0
        The pole position angle of the ellipse in degrees.
        Zero is in the North direction ('g-positive'). Positive clockwise.

    Returns
    -------
    x_model : `float`, array
        Cartesian x-component, in km

    y_model : `float`, array
            Cartesian y-component, in km

    r_model : `float`, array
            Radial distance, in km

    theta : `float` array
            Angular coordinate, in degrees, to return the ellipse points.
    """
    a = equatorial_radius
    b = equatorial_radius - equatorial_radius * oblateness
    phi = position_angle * (np.pi / 180.0)
    ang = theta + phi
    r_model = (a * b) / np.sqrt((a * np.sin(ang)) ** 2 + (b * np.cos(ang)) ** 2)
    x_model = r_model * np.cos(theta) + center_f
    y_model = r_model * np.sin(theta) + center_g
    return x_model, y_model, r_model, theta
