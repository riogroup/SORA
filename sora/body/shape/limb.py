import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
import shapely.geometry as geometry

__all__ = ['Limb']

zero = geometry.Point(0, 0)


class Limb(geometry.LineString):

    @property
    def xy(self):
        return np.array(super(Limb, self).xy)

    def plot(self, center_f=0, center_g=0, scale=1, ax=None, **kwargs):
        """Plots the limb on the tangent plane

        Parameters
        ----------
        center_f : `int`, `float`
            The center of the limb in the x direction. Default=0

        center_g : `int`, `float`
            The center of the limb in the y direction. Default=0

        scale : `int`, `float`
            Scale of the limb relative to the center. Default=1

        ax : `matplotlib.pyplot.Axes`
            The axes where to make the plot. If None, it will use the default axes.

        **kwargs
            All other parameters will be parsed directly by matplotlib.pyplot.

        """
        ax = ax or plt.gca()
        ax.axis('equal')

        xy = self.xy*scale + np.array([[center_f], [center_g]])
        ax.plot(*xy, **kwargs)

    @property
    def maxdist(self):
        """Computes the maximum distance of the limb from the origin"""
        return self.hausdorff_distance(zero)

    def radial_residual_to(self, fg):
        """Calculates radial residuals from points.

        Parameters
        ----------
        fg : `numpy.array`
            Matrix nx2 with the `xy` coordinates of each of the `n` points
            to fit the limb. See example below.

        Returns
        -------
        residuals: `numpy.array`
            Radial distances between limb and given points.

        Examples
        ________

        fg = np.array([[-107.3,   57.8],
                       [ 103.7,   53.2],
                       [ -20.9,  172.4],
                       [   1.9,  171.9]])
        """
        points = geometry.MultiPoint(fg)
        endlines = (fg.T * 1.1 * self.maxdist / np.linalg.norm(fg, axis=-1)).T
        vals = [point.distance(self.intersection(geometry.LineString([zero, endline]))) for point, endline in
                zip(points.geoms, endlines)]
        return np.array(vals)


def limb_radial_residual(limb, fg, center_f=0, center_g=0, scale=1, position_angle=0):
    """

    Parameters
    ----------
    limb : `sora.body.shape.Limb`
        Generic limb to fit.
    fg : `numpy.array`
        Matrix nx2 with the `xy` coordinates of each of the `n` points
        to fit the limb. See example below.
    center_f : `int`, `float`, default=0
        The coordinate in f of the ellipse center.
    center_g : `int`, `float`, default=0
        The coordinate in g of the ellipse center.
    scale : `number`
         Scale factor of the limb
    position_angle : `number`
        The pole position angle of the ellipse in degrees.
        Zero is in the North direction ('g-positive'). Positive clockwise.

    Returns
    -------
    residuals: `numpy.array`
            Radial distances between limb and given points.
    """
    xy = fg - np.array([[center_f], [center_g]]).T
    xy /= scale
    pa = u.Quantity(position_angle, unit='deg')
    rot_mat = np.array([[np.cos(pa), -np.sin(pa)], [np.sin(pa), np.cos(pa)]])
    xy = np.matmul(rot_mat, xy.T).T
    vals = limb.radial_residual_to(xy)
    return np.array(vals)*scale
