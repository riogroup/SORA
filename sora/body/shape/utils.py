import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.coordinates.matrix_utilities import rotation_matrix

__all__ = ['read_obj_file', 'rotated_matrix']


def read_obj_file(filename):
    """Reads a Wavefront OBJ file to get the vertices and faces.

    Parameters
    ----------
    filename : `str`
        Path to the OBJ file.

    Returns
    -------
    vertices : `numpy.array`
        2D array (n, 3) with the "n" vertices of the object

    faces : `numpy.array`
        2D array (m, k) with the "k" vertices that makes each of the "m" faces of the object.

    """
    with open(filename, 'r') as f:
        lines = f.readlines()

    vertices = []
    faces = []
    for line in lines:
        if line.startswith('v '):
            vertices.append(line.strip().split(' ')[1:4])
        elif line.startswith('f '):
            values = line.strip().split(' ')[1:]
            faces.append([i.split('/')[0] for i in values])

    vertices = np.array(vertices, dtype='float')
    faces = np.array(faces, dtype='int32')
    return vertices, faces


def parse_coordinate(coordinate):
    if isinstance(coordinate, str):
        coordinate = SkyCoord(coordinate, unit=(u.deg, u.deg))
    if not isinstance(coordinate, SkyCoord):
        raise ValueError('Coordinate is not in the correct format: (SkyCoord, "+00 00 00.0 +00 00 00.0")')
    return coordinate.spherical


def rotated_matrix(coordinate, sub_observer, pole_position_angle):
    """
    Parameters
    ----------
    coordinate : `astropy.coordinate.CartesianRepresentation`
        Coordinate to be rotated by given orientation
    sub_observer : `astropy.coordinates.SkyCoord`, `str`
        Planetocentric coordinates of the center of the object as seen by the observer.
        It can be an astropy SkyCoord object or a string with the bodycentric longitude
        latitude in degrees. Ex: "30.0 -20.0", or "30 00 00 -20 00 00".
    pole_position_angle : `float`, `int`
        Body's North Pole position angle with respect to direction of the ICRS
        North Pole, i.e. N-E-S-W.

    Returns
    -------
    coordinate : `astropy.coordinate.CartesianRepresentation`
        Coordinate rotated by given orientation
    """
    sub_observer = parse_coordinate(sub_observer)
    pa = u.Quantity(pole_position_angle, unit=u.deg)
    rz = rotation_matrix(-sub_observer.lon, axis='z')
    ry = rotation_matrix(-sub_observer.lat, axis='y')
    rx = rotation_matrix(-pa, axis='x')
    rot_matrix = np.matmul(rx, np.matmul(ry, rz))
    return coordinate.transform(rot_matrix)
