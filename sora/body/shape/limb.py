import matplotlib.pyplot as plt
import numpy as np
from shapely.geometry.polygon import LineString

__all__ = ['Limb']


class Limb(LineString):

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
