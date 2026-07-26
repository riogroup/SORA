import astropy.units as u
import numpy as np
from astropy.coordinates import Angle
from astropy.time import TimeDelta

__all__ = ['Precession']


class Precession:
    """Represent sinusoidal precession terms.

    Parameters
    ----------
    params : `float`, `list`, `numpy.ndarray`
        Parameters used to form the precession equation. It can be a list of
        lists, where each nested list is an independent equation. For example:
        ``params = [1, 2, 3]`` or ``params = [[1, 2], [3, 4]]``. The first
        parameter of each sequence is the amplitude of the sinusoidal function,
        the second parameter is the phase, and the remaining parameters are the
        coefficients of a polynomial time argument.

    func : `str`, default='sin'
        Function applied to the equation. It can be ``'sin'`` for sine or
        ``'cos'`` for cosine.

    multiplier : `str`, default='T'
        Time unit used in the polynomial argument. It can be ``'d'`` for days or
        ``'T'`` for centuries.

    Examples
    --------
    Example with one equation:

    >>> params = [1, 2, 3, 4]
    >>> p = Precession(params, 'sin', 'T')

    The equation to be calculated is:

    >>> val = 1*sin(2 + 3*T + 4*T**2)

    Example with two equations:

    >>> params = [[1, 2, 3, 4, 5], [10, 20, 30, 40, 50]]
    >>> p = Precession(params, 'cos', 'd')

    >>> val = 1*cos(2 + 3*d + 4*d**2 + 5*d**3) + 10*cos(20 + 30*d + 40*d**2 + 50*d**3)

    Notes
    -----
    To inspect the equation, print the `Precession` object.

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
        """Evaluate the precession terms at a time offset.

        Parameters
        ----------
        dt : `float`, `astropy.time.TimeDelta`, `astropy.units.Quantity`
            Time variation from the initial epoch.

        Returns
        -------
        value : `astropy.units.Quantity`
            Total propagated value in degrees.

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
        """Return precession parameters propagated to a time offset.

        Parameters
        ----------
        dt : `float`, `astropy.time.TimeDelta`, `astropy.units.Quantity`
            Time variation from the initial epoch.

        Returns
        -------
        params : `list`
            Input parameters propagated to the requested time.

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
        """Return the Astropy representation of the precession."""
        i, j = self.params.shape
        return "<Precession: {} equations of order {}>".format(i, j - 2)

    def __repr__(self):
        """Return the string representation of the precession."""
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
