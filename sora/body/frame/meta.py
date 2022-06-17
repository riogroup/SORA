import astropy.units as u
import numpy as np
from astropy.coordinates import Angle
from astropy.time import TimeDelta

__all__ = ['Precession']


class Precession:
    """

    Parameters
    ----------
    params: `float`, `list`, `numpy.array`
        List of parameters to form the equation to be calculated.
        It can be a list of list, where each list will be an independent equation.
        For instance: params = [1, 2, 3] or params = [[1, 2],[3, 4]].
        The first parameter of each sequence will be the amplitude of the sinusoidal function.
        The second parameter will be the phase of the sinusoidal function. From the third
        parameter onward, each parameter will be the argument of a linear function.
        See Section "Examples" for more details.

    func: `str`
        Defines the function applied to the equation:
        It can be `sin` for sine and `cos` for cosine.
        Default: `sin`

    multiplier: `str`
        Defines the time unit which the parameters are to be multiplied.
        It can be `d` for day or `T` for centuries.
        Default: `T`

    Examples
    --------

    1)
    >>> params = [1, 2, 3, 4]
    >>> p = Precession(params, 'sin', 'T')

    The equation to be calculated will be

    >>> val = 1*sin(2 + 3*T + 4*T**2)

    2)
    >>> params = [[1, 2, 3, 4, 5], [10, 20, 30, 40, 50]]
    >>> p = Precession(params, 'cos', 'd')

    >>> val = 1*cos(2 + 3*d + 4*d**2 + 5*d**3) + 10*cos(20 + 30*d + 40*d**2 + 50*d**3)

    Notes
    _____
    To chech the real equation just print the Precession object.

    """
    def __init__(self, params=0, func='sin', multiplier='T'):
        if isinstance(params, Precession):
            func = params.func
            multiplier = params.multiplier
            params = params.params
        self.params = np.array(params, ndmin=2)
        if func not in ['sin', 'cos']:
            raise ValueError("'func' must be 'sin' for sine or 'cos' for cosine")
        self.func = func
        if multiplier not in ['d', 'T']:
            raise ValueError("'multiplier' must be 'd' for days or 'T' for centuries.")
        self.multiplier = multiplier
        self.order = len(self.params[0] - 2)

    def compute_at(self, dt):
        """

        Parameters
        ----------
        dt: `float`, `astropy.time.TimeDelta`, `astropy.unit.Quantity`
            Variation in time from the initial epoch.

        Returns
        -------
        : `astropy.unit.Quantity`
            The total value of the parameter propagated to the time required.
            Units: degrees

        """
        dt = TimeDelta(dt)
        func = {'sin': np.sin, 'cos': np.cos}[self.func]
        multi = {'d': u.day, 'T': 100 * u.year}[self.multiplier]
        t = (dt / multi).decompose().value
        tot = 0
        for params in self.params:
            v = np.sum([elem * (t ** i) for i, elem in enumerate(params[1:])]) * u.deg
            tot += params[0] * func(v)
        return tot * u.deg

    def params_at(self, dt):
        """

        Parameters
        ----------
        dt: `float`, `astropy.time.TimeDelta`, `astropy.unit.Quantity`
            Variation in time from the initial epoch.

        Returns
        -------
        : `numpy.array`
            The parameters input propagated to the time required.
            Units: degrees

        """
        dt = TimeDelta(dt)
        multi = {'d': u.day, 'T': 100 * u.year}[self.multiplier]
        t = (dt / multi).decompose().value
        p = []
        for params in self.params:
            v = Angle(np.sum([elem * (t ** i) for i, elem in enumerate(params[1:])]) * u.deg)
            par = [v.wrap_at(360 * u.deg).deg if i == 1 else elem for i, elem in enumerate(params)]
            p.append(par)
        return p

    def _astropy_repr_in_frame(self):
        """Astropy representation of the Precession
        """
        i, j = self.params.shape
        return "<Precession: {} equations of order {}>".format(i, j - 2)

    def __repr__(self):
        """String representation of the Precession Class
        """
        string = []
        m = self.multiplier
        for params in self.params:
            expression = ['{:+f}{}{}'.format(elem, f'*{m}' if i > 0 else "", i if i > 1 else "") for i, elem in
                          enumerate(params[1:])]
            expression = '{}({})'.format(self.func, ''.join(expression)) if len(expression) > 0 else ''
            string.append('{:+f}{}'.format(params[0], expression))
        return '\n'.join(string)

    def __str__(self):
        return self.__repr__()
