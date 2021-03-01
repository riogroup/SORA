import numpy as np

__all__ = ['get_ellipse_points']


def get_ellipse_points(theta, equatorial_radius, oblateness=0.0, center_f=0.0, center_g=0.0,
                       position_angle=0.0):
    """ Get points for the ellipse with the given input parameters

    Parameters:
        theta (float, array): Angular coordinate, in degrees, to return the ellipse points.
        equatorial_radius (float, int): Semi-major axis of the ellipse.
        oblateness (float, int): Oblateness of the ellipse. Default=0.0 (circle)
        center_f (float, int): Coordinate of the ellipse (abscissa). Default=0.0
        center_g (float, int): Coordinate of the ellipse (ordinate). Default=0.0
        position_angle (float, int): The pole position angle of the ellipse in degrees. Default=0
                Zero is in the North direction ('g-positive'). Positive clockwise.
    """
    a = equatorial_radius
    b = equatorial_radius - equatorial_radius * oblateness
    phi = position_angle * (np.pi / 180.0)
    ang = theta + phi
    r_model = (a * b) / np.sqrt((a * np.sin(ang)) ** 2 + (b * np.cos(ang)) ** 2)
    x_model = r_model * np.cos(theta) + center_f
    y_model = r_model * np.sin(theta) + center_g
    return x_model, y_model, r_model, theta
