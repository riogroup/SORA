from functools import lru_cache
from importlib.resources import as_file, files

import numpy as np
import astropy.units as u
import shapely.geometry as geometry
from astropy.coordinates import SkyCoord, CartesianRepresentation
from shapely.ops import unary_union

from sora.config import get_config

from .limb import Limb
from .meta import BaseShape
from .utils import read_obj_file, rotated_matrix

__all__ = ['Shape3D', 'Ellipsoid']


class Shape3D(BaseShape):
    """Represent a 3D shape model loaded from a Wavefront OBJ file.

    Parameters
    ----------
    obj_file : `str`
        Path to the Wavefront OBJ file.
    texture : `str`, optional
        Path to the image that contains the surface texture of the object.
        If `None`, the configured default surface color is used.
    scale : `float`, optional
        Scale factor applied to the shape model vertices. Default is 1.
    right_hand : `bool`, optional
        If `True`, use a right-handed longitude convention. Default is
        `False`.
    """

    def __init__(self, obj_file, texture=None, scale=1, right_hand=False) -> None:
        super(Shape3D, self).__init__()
        vertices, faces = read_obj_file(obj_file)
        self.name = obj_file
        self._vertices = CartesianRepresentation(*vertices.T, unit=u.km)
        self.scale = scale
        if faces.min() == 1:
            faces = faces - 1
        self._faces = faces
        if texture:
            self.texture = texture
        self._right_hand = right_hand

    @property
    def faces(self):
        """`numpy.array` : Vertex indices that define each face of the shape."""
        return self._faces

    @property
    def texture(self):
        """`numpy.array` : RGB colors assigned to each face of the shape."""
        default_color = np.asarray(
            get_config().body_plot.default_surface_color,
            dtype=float,
        )
        texture = np.tile(default_color, (len(self.faces), 1))
        texture = getattr(self, '_texture', texture)
        return texture

    @texture.setter
    def texture(self, value):
        import matplotlib.pyplot as plt
        img = plt.imread(value)
        a, b = img.shape[0:2]
        ym, xm = np.indices(img.shape[0:2])
        yt = (ym + 0.5 - a / 2.0) / (a + 1)
        xt = (xm + 0.5 - b / 2.0) / (b + 1)
        k = np.where(xt < 0.0)
        xt[k] = xt[k] + 1.0

        long_img = -xt * 360 * u.deg + 360 * u.deg
        lat_img = -yt * 180 * u.deg

        img_lats = SkyCoord(ra=-long_img.flatten(), dec=lat_img.flatten())
        faces = SkyCoord(self.vertices[self.faces].mean(axis=-1))

        idx, d2d, d3d = faces.match_to_catalog_sky(img_lats)
        texture = img[ym.flatten()[idx], xm.flatten()[idx]]
        if len(img.shape) == 2:
            texture = np.vstack((texture, texture, texture)).T
        self._texture = np.array(texture) / 256.0

    @property
    def vertices(self):
        """Return the scaled vertices of the shape.

        Returns
        -------
        vertices : `astropy.coordinates.CartesianRepresentation`
            The vertices of the shape.
        """
        return self._vertices * self.scale

    @property
    def scale(self):
        """`float` : Scale factor applied to the shape model vertices."""
        return self._scale

    @scale.setter
    def scale(self, value):
        self._scale = float(value)
        self.get_limb.cache_clear()

    @property
    def volume(self):
        v = self.vertices[self.faces.T]
        ab = v[1] - v[0]
        ac = v[2] - v[0]
        normal = ab.cross(ac)
        return normal.dot(v[0]).sum() / 6

    @property
    def volumetric_equivalent_radius(self):
        return (3 * self.volume / (4 * np.pi)) ** (1. / 3.)

    @property
    def surface_area(self):
        v = self.vertices[self.faces.T]
        ab = v[1] - v[0]
        ac = v[2] - v[0]
        normal = ab.cross(ac)
        return normal.norm().sum() / 2

    @property
    def surface_equivalent_radius(self):
        return np.sqrt(self.surface_area / (4 * np.pi))

    def rotated_vertices(self, sub_observer="00 00 00 +00 00 00", pole_position_angle=0):
        """Return the vertices rotated as viewed by a given observer,
        with 'x' in the direction of the observer.

        Parameters
        ----------
        sub_observer : `astropy.coordinates.SkyCoord`, `str`
            Planetocentric coordinates of the center of the object as seen by the observer.
            It can be an astropy `SkyCoord` object or a string with the
            bodycentric longitude and latitude in degrees. For example:
            ``"30.0 -20.0"`` or ``"30 00 00 -20 00 00"``.

        pole_position_angle : `float`, `int`
            Body's North Pole position angle with respect to direction of the ICRS
            North Pole, i.e. N-E-S-W.

        Returns
        -------
        rotated_vertices : `astropy.coordinates.CartesianRepresentation`
            An astropy object with the rotated vertices as viewed by the observer.
            The 'x' coordinate is in the direction of the observer, the 'z' coordinate
            in the direction of the projected ICRS North Pole, and the 'y' coordinate
            complementing the right-hand rule.
        """
        rotated_vertices = rotated_matrix(coordinate=self.vertices, sub_observer=sub_observer,
                                          pole_position_angle=pole_position_angle, right_hand=self._right_hand)
        return rotated_vertices

    @lru_cache(maxsize=128)
    def get_limb(self, sub_observer="00 00 00 +00 00 00", pole_position_angle=0, center_f=0, center_g=0, **kwargs):
        """Return the projected limb as viewed by a given observer,
        with 'x' in the direction of the observer.

        Parameters
        ----------
        sub_observer : `astropy.coordinates.SkyCoord`, `str`
            Planetocentric coordinates of the center of the object as seen by the observer.
            It can be an astropy `SkyCoord` object or a string with the
            bodycentric longitude and latitude in degrees. For example:
            ``"30.0 -20.0"`` or ``"30 00 00 -20 00 00"``.

        pole_position_angle : `float`, `int`
            Body's North Pole position angle with respect to direction of the ICRS
            North Pole, i.e. N-E-S-W.

        center_f : `int`, `float`
            Offset of the center of the body in the East direction, in km.

        center_g : `int`, `float`
            Offset of the center of the body in the North direction, in km.

        **kwargs
            Any other keyword argument is discarded.

        Returns
        -------
        limb : `sora.body.shape.Limb`
            The limb corresponding to the projected view of the 3D shape.
        """
        vertices = self.rotated_vertices(sub_observer=sub_observer, pole_position_angle=pole_position_angle)
        faced_vertices = vertices[self.faces.T]

        ab = faced_vertices[0] - faced_vertices[1]
        ac = faced_vertices[0] - faced_vertices[-1]
        normal_obs = ab.cross(ac)
        observable = normal_obs.x >= 0

        x = - faced_vertices.T[observable].y.value + center_f
        y = faced_vertices.T[observable].z.value + center_g
        triangles = [geometry.Polygon(np.transpose(triangle)) for triangle in zip(x, y)]
        
        valid_triangles = [poly.buffer(0) if not poly.is_valid else poly for poly in triangles]
        
        pol = unary_union(valid_triangles)
        limb = pol.boundary
        if isinstance(limb, geometry.multilinestring.MultiLineString):
            limb = limb.geoms[0]
        return Limb(limb)

    def plot(self, sub_observer="00 00 00 +00 00 00", sub_solar=None,
             pole_position_angle=0, center_f=0, center_g=0, scale=1, ax=None,
             plot_pole=None, north_pole_color=None, south_pole_color=None,
             pole_length_scale=None, surface_alpha=None, **kwargs):
        """Plot the projected 3D shape on the tangent plane.

        Parameters
        ----------
        sub_observer : `astropy.coordinates.SkyCoord`, `str`
            Planetocentric coordinates of the center of the object in the direction of the observer.
            It can be an astropy `SkyCoord` object or a string with the
            bodycentric longitude and latitude in degrees. For example:
            ``"30.0 -20.0"`` or ``"30 00 00 -20 00 00"``.

        sub_solar : `astropy.coordinates.SkyCoord`, `str`
            Planetocentric coordinates of the center of the object in the direction of the Sun.
            It can be an astropy `SkyCoord` object or a string with the
            bodycentric longitude and latitude in degrees. For example:
            ``"30.0 -20.0"`` or ``"30 00 00 -20 00 00"``.

        pole_position_angle : `float`, `int`
            Body's North Pole position angle with respect to direction of the ICRS
            North Pole, i.e. N-E-S-W.

        center_f : `int`, `float`
            Offset of the center of the body in the East direction, in km.

        center_g : `int`, `float`
            Offset of the center of the body in the North direction, in km.

        scale : `float`
            Factor applied to the plotted shape vertices. Default is 1.

        ax : `matplotlib.pyplot.Axes`
            The axes where to make the plot. If None, it will use the default axes.

        plot_pole : `bool`, optional
            If True, the direction of the pole is plotted. When omitted, uses
            the configured value.

        north_pole_color : `str`, optional
            Color used to plot the North pole. When omitted, uses the
            configured value.

        south_pole_color : `str`, optional
            Color used to plot the South pole. When omitted, uses the
            configured value.

        pole_length_scale : `float`, optional
            Scale applied to the plotted pole length. When omitted, uses the
            configured value.

        surface_alpha : `float`, optional
            Opacity of the plotted shape surface, between 0 and 1. When
            omitted, uses the configured value.

        **kwargs
            Any other keyword argument is discarded.
        """
        import matplotlib.pyplot as plt

        plot_config = get_config().body_plot
        if plot_pole is None:
            plot_pole = plot_config.show_pole
        if north_pole_color is None:
            north_pole_color = plot_config.north_pole_color
        if south_pole_color is None:
            south_pole_color = plot_config.south_pole_color
        if pole_length_scale is None:
            pole_length_scale = plot_config.pole_length_scale
        if surface_alpha is None:
            surface_alpha = plot_config.surface_alpha

        ax = ax or plt.gca()
        ax.axis('equal')

        params = {'sub_observer': sub_observer, 'pole_position_angle': pole_position_angle}
        vertices_observer = self.rotated_vertices(**params)
        cart_obs = vertices_observer[self.faces.T]
        ab = cart_obs[0] - cart_obs[1]
        ac = cart_obs[0] - cart_obs[-1]
        normal_obs = ab.cross(ac)
        observable = normal_obs.x >= 0

        if not sub_solar:
            sub_solar = sub_observer
        vertices_sun = self.rotated_vertices(sub_observer=sub_solar)

        cart_sun = vertices_sun[self.faces.T]
        ab = cart_sun[0] - cart_sun[1]
        ac = cart_sun[0] - cart_sun[2]
        normal_sun = ab.cross(ac)
        normal_sun = normal_sun / normal_sun.norm()

        shade = normal_sun.x
        shade[shade < 0] = 0
        alpha = np.full(len(shade), surface_alpha)
        color = np.vstack((self.texture.T * shade, alpha)).T.value

        # Define pole for plotting
        vert = self.vertices
        maxd = vert.norm().max()
        maxz = vert.z.argmax()
        minz = vert.z.argmin()
        npole = CartesianRepresentation(
            [0, 0] * u.km,
            [0, 0] * u.km,
            [vert[maxz].norm(), maxd * pole_length_scale],
        ) * scale
        spole = CartesianRepresentation(
            [0, 0] * u.km,
            [0, 0] * u.km,
            [-vert[minz].norm(), -maxd * pole_length_scale],
        ) * scale
        params['right_hand'] = self._right_hand
        rnpole = rotated_matrix(coordinate=npole, **params)
        rspole = rotated_matrix(coordinate=spole, **params)

        poles = {'north': rnpole, 'south': rspole}
        pcolor = {
            'north': north_pole_color,
            'south': south_pole_color,
        }
        forepole = 'north' if rnpole[0].x > rspole[0].x else 'south'
        backpole = 'south' if forepole == 'north' else 'north'

        # Plots the pole that is in the background
        if plot_pole:
            ax.plot(-poles[backpole].y.value + center_f, poles[backpole].z.value + center_g, color=pcolor[backpole], zorder=1)

        for i, pol in enumerate(scale*cart_obs.T):
            if not observable[i]:
                continue
            ax.fill(-pol.y.value + center_f, pol.z.value + center_g, color=color[i], zorder=1)

        # Plots the pole that is in the foreground
        if plot_pole:
            ax.plot(-poles[forepole].y.value + center_f, poles[forepole].z.value + center_g, color=pcolor[forepole], zorder=1)


class Ellipsoid(Shape3D):
    """Represent a triaxial ellipsoid as a 3D shape model.

    Parameters
    ----------
    a : `float`
        Semi-axis along the x direction, in km.
    b : `float`, optional
        Semi-axis along the y direction, in km. If `None`, it is set to `a`.
    c : `float`, optional
        Semi-axis along the z direction, in km. If `None`, it is set to `b`.
    texture : `str`, optional
        Path to the image that contains the surface texture of the object.
        If `None`, a gray color is used.
    right_hand : `bool`, optional
        If `True`, use a right-handed longitude convention. Default is
        `False`.
    """

    def __init__(self, a, b=None, c=None, texture=None, right_hand=False):
        sphere = files('sora').joinpath('data/sphere.obj')
        with as_file(sphere) as obj_file:
            super().__init__(obj_file=str(obj_file), texture=texture, right_hand=right_hand)
        self.a = a
        self.b = b or self.a
        self.c = c or self.b
        self.name = f'{self.a} x {self.b} x {self.c}'
        v = self.vertices
        norm = v/np.sqrt(v.dot(v)).value
        self._vertices = CartesianRepresentation(norm.x*self.a, norm.y*self.b, norm.z*self.c)

    @lru_cache(maxsize=128)
    def get_limb(self, sub_observer="00 00 00 +00 00 00", pole_position_angle=0, center_f=0, center_g=0):
        # TODO(Compute the limb from equation to avoid unnecessary use of the 3D shape)
        return super(Ellipsoid, self).get_limb(sub_observer=sub_observer, pole_position_angle=pole_position_angle,
                                               center_f=center_f, center_g=center_g)

    get_limb.__doc__ = Shape3D.get_limb.__doc__
