import numpy as np

__all__ = ['ChiSquare']


class ChiSquare:
    """Stores the arrays for all inputs and given chi-square.

    Parameters
    ----------
    chi2 : `array`
        Array with all the chi-square values.

    npts : `int`
        Number of points used in the fit.

    **kwargs
        Any other given input must be an array with the same size as chi2. The
        keyword `name` will be associated as the variable `name` of the given data.

    Example
    -------
    >>> chisquare = ChiSquare(chi2, immersion=t1, emersion=t2)

    ``t1`` and ``t2`` must be an array with the same size as chi2.

    The data can be accessed as

    >>> chisquare.data['immersion']

    """

    def __init__(self, chi2, npts, **kwargs):

        self._names = ['chi2']
        self.data = {'chi2': chi2}
        data_size = len(chi2)
        self.npts = npts
        nparam = 0
        for item in kwargs.keys():
            if kwargs[item].var() != 0:
                nparam += 1
            if len(kwargs[item]) != data_size:
                raise ValueError('{} size must have the same size as given chi2'.format(item))
            self._names.append(item)
            self.data[item] = kwargs[item]
        self.nparam = nparam

    def get_nsigma(self, sigma=1, key=None):
        """Determines the interval of the chi-square within the nth sigma.

        Parameters
        ----------
        sigma : `float`, `int`
            Value of sigma to calculate.

        key : `str`, default=None
            keyword the user desire to obtain results.

        Returns
        -------
        dict
            Dictionary with the average n-sigma value and bondaries.

        Note
        ----
        If a key value is given the mean value within the n-sigma and the error
        bar within the n-sigma are returned.

        If no key is given, a dictionary with: the minimum chi-square, the sigma
        required, the number of points where :math:`chi2 < (chi2_min + sigma^2)`,
        and the mean values and errors for all keys is returned.
        """
        values = np.where(self.data['chi2'] < self.data['chi2'].min() + sigma ** 2)[0]
        output = {'chi2_min': self.data['chi2'].min(), 'sigma': sigma, 'n_points': len(values)}
        for name in self._names[1:]:
            vmax = self.data[name][values].max()
            vmin = self.data[name][values].min()
            output[name] = [(vmax + vmin) / 2.0, (vmax - vmin) / 2.0]
        if key is not None:
            if key not in self._names[1:]:
                raise ValueError('{} is not one of the available keys. Please choose one of {}'
                                 .format(key, self._names[1:]))
            return output[key]
        return output

    def plot_chi2(self, key=None, ax=None):
        """Plots an ellipse using the input parameters.

        Parameters
        ----------
        key : `str`
            Key (parameter) for which to plot the chi squares. If no key  is
            given, it will plot for all parameters.

        ax : `maplotlib.pyplot.Axes`
            Matplotlib plot axes. If none is given, it uses the matplotlib pool
            to identify.
        """
        import matplotlib.pyplot as plt

        sigma_1 = self.get_nsigma(sigma=1)
        sigma_3 = self.get_nsigma(sigma=3)
        if key is not None and (key not in self._names[1:]):
            raise ValueError('{} is not one of the available keys. Please choose one of {}'
                             .format(key, self.data.keys()[1:]))
        k = np.where(self.data['chi2'] < sigma_1['chi2_min'] + 11)
        for name in self._names[1:]:
            if (key is not None) and (key != name):
                continue
            if key is None or not ax:
                ax = plt.gca()
            ax.plot(self.data[name][k], self.data['chi2'][k], 'k.')
            ax.set_ylim(sigma_1['chi2_min'] - 1, sigma_1['chi2_min'] + 10)
            delta = sigma_3[name][1]
            if delta == 0.0:
                delta = 1.0
            ax.set_xlim(sigma_3[name][0] - 3 * delta, sigma_3[name][0] + 3 * delta)
            ax.axhline(sigma_1['chi2_min'] + 1, linestyle='--', color='red')
            ax.axhline(sigma_3['chi2_min'] + 9, linestyle='--', color='red')
            ax.set_xlabel(name, fontsize=20)
            ax.set_ylabel('Chi2', fontsize=20)
            if key is None:
                plt.show()

    def to_file(self, namefile):
        """Saves the data to a file.

        Parameters
        ----------
        namefile : `str`
            Filename to save the data.
        """
        data = np.vstack(([self.data[i] for i in self._names]))
        np.savetxt(namefile, data.T, fmt='%11.5f')
        f = open(namefile + '.label', 'w')
        for i, name in enumerate(self._names):
            f.write('Column {}: {}\n'.format(i + 1, name))
        f.close()

    def get_values(self, sigma=0.0, key=None):
        """Returns all values where the chi-square is within the nth sigma.

        Parameters
        ----------
        sigma : `float`, `int`
            Value of sigma to cut values.

        key : `str`
            Keyword (parameter) that the user desires to obtain results of.

        Returns
        -------
        list or dict : `list`, `dict`
            List or dictionary with chi-square values within the nth sigma, the
            average n-sigma value, and bondaries.

        Note
        ----
        If a `key` is given, it returns a `list` with all the values that are within
        the chosen n-sigma.

        If no `key` is given, it returns a dictionary with the list with all the
        values that are within the n-sigma for all keys.

        If ``sigma=0``, it returns the parameters for the minimum chi-square
        instead of a list.
        """
        values = {}
        if sigma == 0.0:
            k = np.argsort(self.data['chi2'])[0]
        else:
            k = np.where(self.data['chi2'] < self.data['chi2'].min() + sigma ** 2)[0]
        for name in self._names[1:]:
            values[name] = self.data[name][k]
        values = values.get(key, values)
        return values

    def __len__(self):
        return len(self.data['chi2'])

    def __add__(self, other):
        if not isinstance(other, ChiSquare):
            raise TypeError(
                f"unsupported operand type(s) for +: '{self.__class__.__name__}' and '{other.__class.__name}'")
        if self._names != other._names:
            raise ValueError(f"ChiSquare objects does not have the same keys: '{self._names}' and '{other._names}'")
        if self.npts != other.npts:
            raise ValueError(
                f"The number of fitted points are different between the objects: '{self.npts}' and '{other.npts}'")
        params = {key: np.hstack((self.data[key], other.data[key])) for key in self._names}
        chi2 = params.pop('chi2')
        return ChiSquare(chi2=chi2, npts=int((self.npts + other.npts) / 2), **params)

    def __str__(self):
        """ String representation of the ChiSquare Object
        """
        sigma_1 = self.get_nsigma(sigma=1)
        sigma_3 = self.get_nsigma(sigma=3)
        output = ('Minimum chi-square: {:.3f}\n'
                  'Number of fitted points: {}\n'
                  'Number of fitted parameters: {}\n'
                  'Minimum chi-square per degree of freedom: {:.3f}\n'.format(
                      sigma_1['chi2_min'], self.npts, self.nparam,
                      sigma_1['chi2_min']/(self.npts-self.nparam))
                  )
        for name in self._names[1:]:
            output += ('\n{}:\n'
                       '    1-sigma: {:.3f} +/- {:.3f}\n'
                       '    3-sigma: {:.3f} +/- {:.3f}\n'.format(
                           name, *sigma_1[name], *sigma_3[name])
                       )
        return output
