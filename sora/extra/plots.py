import astropy.units as u
import numpy as np

__all__ = ['draw_ellipse']


def draw_ellipse(equatorial_radius, oblateness=0.0, center_f=0.0, center_g=0.0,
                 position_angle=0.0, center_dot=False, ax=None, **kwargs):
    """Plots an ellipse with the given input parameters.

    Parameters
    ----------
    equatorial_radius : `float`, `int`, default=0
        Semi-major axis of the ellipse.

    oblateness : `float`, `int`, default=0
        Oblateness of the ellipse.

    center_f : `float`, `int`, default=0
        Coordinate of the ellipse (abscissa).

    center_g : `float`, `int`, default=0
        Coordinate of the ellipse (ordinate).

    center_dot : `bool`, default=False
        If True, plots a dot at the center of the ellipse.

    position_angle : `float`, `int`, default= 0
        Pole position angle (Default=0.0).

    ax : `maptlotlib.pyplot.Axes`
        Axis where to plot ellipse.

    **kwargs
        All other parameters. They will be parsed directly by matplotlib.
    """
    import matplotlib.pyplot as plt

    equatorial_radius = np.array(equatorial_radius, ndmin=1)
    oblateness = np.array(oblateness, ndmin=1)
    center_f = np.array(center_f, ndmin=1)
    center_g = np.array(center_g, ndmin=1)
    position_angle = np.array(position_angle, ndmin=1)

    theta = np.linspace(-np.pi, np.pi, 1800)

    ax = ax or plt.gca()

    if len(equatorial_radius) == 1:
        if 'color' not in kwargs:
            kwargs['color'] = 'black'
        if 'lw' not in kwargs:
            kwargs['lw'] = 2
    else:
        if 'color' not in kwargs:
            kwargs['color'] = 'gray'
        if 'lw' not in kwargs:
            kwargs['lw'] = 0.1
        if 'alpha' not in kwargs:
            kwargs['alpha'] = 0.1
        if 'zorder' not in kwargs:
            kwargs['zorder'] = 0.5
    for i in np.arange(len(equatorial_radius)):
        circle_x = equatorial_radius[i] * np.cos(theta)
        circle_y = equatorial_radius[i] * (1.0 - oblateness[i]) * np.sin(theta)
        pos_ang = position_angle[i] * u.deg
        ax.plot(+circle_x * np.cos(pos_ang) + circle_y * np.sin(pos_ang) + center_f[i],
                -circle_x * np.sin(pos_ang) + circle_y * np.cos(pos_ang) + center_g[i],
                **kwargs)
    if center_dot:
        kwargs.pop('lw')
        plt.plot(center_f, center_g, '.', **kwargs)
    plt.axis('equal')
