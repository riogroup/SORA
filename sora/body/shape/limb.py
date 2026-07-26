import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
import shapely.geometry as geometry

__all__ = ['Limb']

zero = geometry.Point(0, 0)


class Limb:
    __doc__ = geometry.LineString.__doc__

    def __init__(self, *args, **kwargs):
        self._contour = geometry.LineString(*args, **kwargs)
        self.shadow = geometry.Polygon(self._contour)

    @property
    def xy(self):
        """`numpy.array` : Coordinates of the limb contour."""
        return np.array(self._contour.xy)

    @property
    def shadow_area(self):
        return self.shadow.area * u.km**2

    @property
    def shadow_equivalent_radius(self):
        return np.sqrt(self.shadow_area/np.pi)

    def plot(self, center_f=0, center_g=0, scale=1, ax=None, **kwargs):
        """Plot the limb on the tangent plane.

        Parameters
        ----------
        center_f : `int`, `float`
            The center of the limb in the f direction. Default is 0.

        center_g : `int`, `float`
            The center of the limb in the g direction. Default is 0.

        scale : `int`, `float`
            Scale of the limb relative to the center. Default is 1.

        ax : `matplotlib.pyplot.Axes`
            The axes where to make the plot. If None, it will use the default axes.

        **kwargs
            All other parameters are passed directly to `matplotlib.pyplot`.

        """
        ax = ax or plt.gca()
        ax.axis('equal')

        xy = self.xy*scale + np.array([[center_f], [center_g]])
        ax.plot(*xy, **kwargs)

    @property
    def maxdist(self):
        """`float` : Maximum distance of the limb from the origin."""
        return self._contour.hausdorff_distance(zero)

    def radial_residual_to(self, fg):
        """Calculate radial residuals from points.

        Parameters
        ----------
        fg : `numpy.array`
            Matrix (n, 2) with the f and g coordinates of each of the `n`
            points to compare with the limb.

        Returns
        -------
        residuals : `numpy.array`
            Radial distances between limb and given points.

        Examples
        --------

        >>> fg = np.array([[-107.3, 57.8],
        ...                [103.7, 53.2],
        ...                [-20.9, 172.4],
        ...                [1.9, 171.9]])
        """
        points = geometry.MultiPoint(fg)
        endlines = (fg.T * 1.1 * self.maxdist / np.linalg.norm(fg, axis=-1)).T
        vals = [point.distance(self._contour.intersection(geometry.LineString([zero, endline]))) for point, endline in
                zip(points.geoms, endlines)]
        return np.array(vals)

    def __str__(self):
        return self._contour.__str__()

    def __repr__(self):
        return self._contour.__repr__()


def limb_radial_residual(limb, fg, center_f=0, center_g=0, scale=1, position_angle=0):
    """Calculate radial residuals after applying limb fit parameters.

    Parameters
    ----------
    limb : `sora.body.shape.Limb`
        Generic limb to fit.

    fg : `numpy.array`
        Matrix (n, 2) with the f and g coordinates of each of the `n`
        points to fit the limb.

    center_f : `int`, `float`, default=0
        The coordinate in f of the limb center.

    center_g : `int`, `float`, default=0
        The coordinate in g of the limb center.

    scale : `number`
        Scale factor of the limb.

    position_angle : `number`
        The pole position angle of the limb in degrees.
        Zero is in the North direction ('g-positive'). Positive clockwise.

    Returns
    -------
    residuals : `numpy.array`
        Radial distances between limb and given points.
    """
    xy = fg - np.array([[center_f], [center_g]]).T
    xy /= scale
    pa = u.Quantity(position_angle, unit='deg')
    rot_mat = np.array([[np.cos(pa), -np.sin(pa)], [np.sin(pa), np.cos(pa)]])
    xy = np.matmul(rot_mat, xy.T).T
    vals = limb.radial_residual_to(xy)
    return np.array(vals)*scale
