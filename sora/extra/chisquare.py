import numpy as np

__all__ = ['ChiSquare']


class ChiSquare:
    """Store chi-square values and associated fitted parameters.

    Parameters
    ----------
    chi2 : `numpy.array`
        Array with all the chi-square values.

    npts : `int`
        Number of points used in the fit.

    **kwargs
        Additional inputs. Each value must be an array with the same size as
        `chi2`. Each keyword is used as the parameter name in `data`.

    Examples
    --------
    >>> chisquare = ChiSquare(chi2, npts=4, immersion=t1, emersion=t2)

    ``t1`` and ``t2`` must be arrays with the same size as ``chi2``.

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
        """Return intervals for values within the requested sigma level.

        Parameters
        ----------
        sigma : `float`, `int`
            Sigma level used to select values with
            ``chi2 < chi2_min + sigma**2``.

        key : `str`, default=None
            Parameter name for which to return the interval. If `None`,
            intervals for all parameters are returned.

        Returns
        -------
        result : `dict`, `list`
            If `key` is `None`, a dictionary with the minimum chi-square,
            sigma level, number of selected points, and intervals for all
            parameters. If `key` is given, a list with the central value and
            half-width for that parameter.

        Raises
        ------
        ValueError
            If `key` is not one of the available parameter names.

        Notes
        -----
        Parameter intervals are computed from the minimum and maximum values
        selected within the requested sigma level.
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
        """Plot chi-square values as a function of fitted parameters.

        Parameters
        ----------
        key : `str`, optional
            Parameter name for which to plot the chi-square values. If `None`,
            a plot is generated for each parameter.

        ax : `matplotlib.pyplot.Axes`, optional
            Matplotlib axes where the plot is drawn. If `None`, the current
            axes are used.

        Raises
        ------
        ValueError
            If `key` is not one of the available parameter names.
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
        """Save chi-square data and column labels to files.

        Parameters
        ----------
        namefile : `str`
            Filename used to save the data. Column labels are written to
            ``namefile + '.label'``.
        """
        data = np.vstack(([self.data[i] for i in self._names]))
        np.savetxt(namefile, data.T, fmt='%11.5f')
        f = open(namefile + '.label', 'w')
        for i, name in enumerate(self._names):
            f.write('Column {}: {}\n'.format(i + 1, name))
        f.close()

    def get_values(self, sigma=0.0, key=None):
        """Return values where chi-square is within the requested sigma level.

        Parameters
        ----------
        sigma : `float`, `int`
            Sigma level used to select values with
            ``chi2 < chi2_min + sigma**2``. If 0, only values at the minimum
            chi-square are returned.

        key : `str`, optional
            Parameter name for which to return values. If `None`, values for
            all parameters are returned.

        Returns
        -------
        values : `dict`, `numpy.array`
            If `key` is `None`, a dictionary with selected values for all
            parameters. If `key` is given, the selected values for that
            parameter.

        Notes
        -----
        If ``sigma=0``, this returns the parameter values at the minimum
        chi-square instead of arrays selected by a sigma threshold.
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

    def to_log(self, namefile):
        """Save the chi-square summary log to a file.

        Parameters
        ----------
        namefile : `str`
            Filename to save the log.
        """
        f = open(namefile, 'w')
        f.write(self.__str__())
        f.close()

    def __len__(self):
        """Return the number of chi-square values stored."""
        return len(self.data['chi2'])

    def __add__(self, other):
        """Concatenate two compatible `ChiSquare` objects.

        Parameters
        ----------
        other : `ChiSquare`
            Object to concatenate with this one.

        Returns
        -------
        chisquare : `ChiSquare`
            New object containing the concatenated chi-square and parameter
            arrays.

        Raises
        ------
        TypeError
            If `other` is not a `ChiSquare` object.
        ValueError
            If the two objects have different parameter names or different
            numbers of fitted points.
        """
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
        """Return a printable summary of the chi-square fit."""
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
