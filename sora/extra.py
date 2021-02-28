import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u
from sora.config.visuals import progressbar


def draw_ellipse(equatorial_radius, oblateness=0.0, center_f=0.0, center_g=0.0,
                 position_angle=0.0, center_dot=False, ax=None, **kwargs):
    """ Plots an ellipse with the given input parameters

    Parameters:
        radius (float, int): Semi-major axis of the ellipse.
        oblateness (float, int): Oblateness of the ellipse. Default=0.0 (circle)
        center_x (float, int): Coordinate of the ellipse (abscissa). Default=0.0
        center_y (float, int): Coordinate of the ellipse (ordinate). Default=0.0
        center_dot (bool): If True, plots a dot at the center of the ellipse. Default=False
        position_angle (float, int): Pole position angle. Default=0.0
        ax (maptlotlib.Axes): Axis where to plot ellipse
        **kwargs: all other parameters will be parsed directly to matplotlib
    """
    equatorial_radius = np.array(equatorial_radius, ndmin=1)
    oblateness = np.array(oblateness, ndmin=1)
    center_f = np.array(center_f, ndmin=1)
    center_g = np.array(center_g, ndmin=1)
    position_angle = np.array(position_angle, ndmin=1)

    theta = np.linspace(-np.pi, np.pi, 1800)
    size_vals, size_theta = np.indices((len(equatorial_radius), len(theta)))

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
        circle_x = equatorial_radius[i]*np.cos(theta)
        circle_y = equatorial_radius[i]*(1.0-oblateness[i])*np.sin(theta)
        pos_ang = position_angle[i]*u.deg
        ax.plot(+circle_x*np.cos(pos_ang) + circle_y*np.sin(pos_ang) + center_f[i],
                -circle_x*np.sin(pos_ang) + circle_y*np.cos(pos_ang) + center_g[i],
                **kwargs)
    if center_dot:
        kwargs.pop('lw')
        plt.plot(center_f, center_g, '.', **kwargs)
    plt.axis('equal')

def get_ellipse_points(theta, equatorial_radius, oblateness=0.0, center_f=0.0, center_g=0.0,
                 position_angle=0.0):
    """ Get points for the ellipse with the given input parameters

    Parameters:
        radius (float, int): Semi-major axis of the ellipse.
        oblateness (float, int): Oblateness of the ellipse. Default=0.0 (circle)
        center_x (float, int): Coordinate of the ellipse (abscissa). Default=0.0
        center_y (float, int): Coordinate of the ellipse (ordinate). Default=0.0
        position_angle (float, int): The pole position angle of the ellipse in degrees. Default=0
                Zero is in the North direction ('g-positive'). Positive clockwise.
    """
    a = equatorial_radius
    b = equatorial_radius - equatorial_radius*oblateness
    phi = position_angle*(np.pi/180.0)
    ang = theta+phi
    r_model = (a*b)/np.sqrt((a*np.sin(ang))**2 + (b*np.cos(ang))**2)
    x_model = r_model*np.cos(theta) + center_f
    y_model = r_model*np.sin(theta) + center_g
    return x_model, y_model, r_model, theta

def filter_negative_chord(chord,chisquare,step=1,sigma=0):
    """ Get points for the ellipse with the given input parameters

    Parameters:
        chord (Chord): Chord object, must be associated to an Occultation to work.
        chisquare (ChiSquare): Resulted ChiSquare object of fit_ellipse.
        sigma (int, float): Unceartity of the expected ellipse, in km.
        step (number, 'exposure'): If a number, it corresponds to the step, in seconds, for each point of the chord path.
            The step can also be equal to 'exposure'. In this case, the chord path will consider the lightcurve individual
            times and exptime.
    """
    keep = []
    if step == 'exposure':
        try:
            step = np.min([chord.lightcurve.exptime/10., step])
        except:
            raise ValueError('Chord.lightcurve does not have "exptime"')
        time_all = np.arange(chord.lightcurve.time.min(), chord.lightcurve.time.max(), step)
        time_exposure = np.array([])
        for i in range(len(chord.lightcurve.time)):
            event_model = (time_all > chord.lightcurve.time[i]-chord.lightcurve.exptime/2.) & (time_all < chord.lightcurve.time[i]+chord.lightcurve.exptime/2.)
            time_exposure = np.append(time_exposure,time_all[event_model])
            f_all, g_all = chord.get_fg(time=time_exposure*u.s + chord.lightcurve.tref)
    else:
        f_all, g_all = chord.path(segment='full', step=step)
    for i in progressbar(range(len(chisquare.data['chi2'])),'Filter chord: {}'.format(chord.name)):
        df_all = (f_all - chisquare.data['center_f'][i])
        dg_all = (g_all - chisquare.data['center_g'][i])

        r_all = np.sqrt(df_all**2 + dg_all**2)
        cut = r_all < 1.5*chisquare.data['equatorial_radius'][i]
        df_path = df_all[cut]
        dg_path = dg_all[cut]
        r_path = r_all[cut]
        theta_path = np.arctan2(dg_path, df_path)

        r_ellipse = get_ellipse_points(theta_path,
                                       equatorial_radius=chisquare.data['equatorial_radius'][i],
                                       oblateness=chisquare.data['oblateness'][i],
                                       center_f=chisquare.data['center_f'][i],
                                       center_g=chisquare.data['center_g'][i],
                                       position_angle=chisquare.data['position_angle'][i])[2]
        keep.append(np.all(r_path - r_ellipse + sigma > 0))

    filtered_chisquare = ChiSquare(chisquare.data['chi2'][keep],chisquare.npts,
                                   center_f= chisquare.data['center_f'][keep],
                                   center_g= chisquare.data['center_g'][keep],
                                   equatorial_radius= chisquare.data['equatorial_radius'][keep],
                                   oblateness= chisquare.data['oblateness'][keep],
                                   position_angle= chisquare.data['position_angle'][keep])
    return filtered_chisquare

class ChiSquare():
    def __init__(self, chi2, npts, **kwargs):
        """ Stores the arrays for all inputs and given chi-square.

        Parameters:
            chi2 (array): Array with all the chi-square values
            npts (int): Number of points used in the fit
            **kwargs: any other given input must be an array with the same size as chi2.
                the keyword name will be associated as the variable name of the given data

        Example:

        chisquare = ChiSquare(chi2, immersion=t1, emersion=t2)
            t1 and t2 must be an array with the same size as chi2.
        the data can be accessed as:
            chisquare.data['immersion']
        """
        self.__names = ['chi2']
        self.data = {'chi2': chi2}
        data_size = len(chi2)
        self.npts = npts
        nparam = 0
        for item in kwargs.keys():
            if (kwargs[item].var() != 0):
                nparam += 1
            if len(kwargs[item]) != data_size:
                raise ValueError('{} size must have the same size as given chi2')
            self.__names.append(item)
            self.data[item] = kwargs[item]
        self.nparam = nparam

    def get_nsigma(self, sigma=1, key=None):
        """ Determines the interval of the chi-square within the n-th sigma

        Parameters:
            sigma (float, int): Value of sigma to calculate.
            key (str): keyword the user desire to obtain results. Default=None

        Return:
            - if a key is given, it returns two values: the mean value within the n-sigma
                and the error bar within the n-sigma.
            - if no key is given, it returns a dictionary with the minimum chi-square,
                the sigma required, the number of points where chi2 < chi2_min + sigma^2,
                and the mean values and errors for all keys.
        """
        values = np.where(self.data['chi2'] < self.data['chi2'].min() + sigma**2)[0]
        output = {'chi2_min': self.data['chi2'].min()}
        output['sigma'] = sigma
        output['n_points'] = len(values)
        for name in self.__names[1:]:
            vmax = self.data[name][values].max()
            vmin = self.data[name][values].min()
            output[name] = [(vmax+vmin)/2.0, (vmax-vmin)/2.0]
        if key is not None:
            if key not in self.__names[1:]:
                raise ValueError('{} is not one of the available keys. Please choose one of {}'
                                 .format(key, self.__names[1:]))
            return output[key]
        return output

    def plot_chi2(self, key=None, ax=None):
        """ Plots an ellipse using the given input parameters

        Parameters:
            key (str): Key to plot chi square.
                if no key  is given, it will plot for all the keywords.
            ax (maplotlib.Axes): A matplotlib axes to plot,
                if none is given, it uses the matplotlib pool to identify.
        """
        sigma_1 = self.get_nsigma(sigma=1)
        sigma_3 = self.get_nsigma(sigma=3)
        if key is not None and (key not in self.__names[1:]):
            raise ValueError('{} is not one of the available keys. Please choose one of {}'
                             .format(key, self.data.keys()[1:]))
        k = np.where(self.data['chi2'] < sigma_1['chi2_min']+11)
        for name in self.__names[1:]:
            if (key is not None) and (key != name):
                continue
            if key is None or not ax:
                ax = plt.gca()
            ax.plot(self.data[name][k], self.data['chi2'][k], 'k.')
            ax.set_ylim(sigma_1['chi2_min']-1, sigma_1['chi2_min']+10)
            delta = sigma_3[name][1]
            if delta == 0.0:
                delta = 1.0
            ax.set_xlim(sigma_3[name][0]-3*delta, sigma_3[name][0]+3*delta)
            ax.axhline(sigma_1['chi2_min']+1, linestyle='--', color='red')
            ax.axhline(sigma_3['chi2_min']+9, linestyle='--', color='red')
            ax.set_xlabel(name, fontsize=20)
            ax.set_ylabel('Chi2', fontsize=20)
            if key is None:
                plt.show()

    def to_file(self, namefile):
        """ Saves the data to a file

        Parameters:
            namefile (str): Filename to save the data
        """
        data = np.vstack(([self.data[i] for i in self.__names]))
        np.savetxt(namefile, data.T, fmt='%11.5f')
        f = open(namefile+'.label', 'w')
        for i, name in enumerate(self.__names):
            f.write('Column {}: {}\n'.format(i+1, name))
        f.close()

    def get_values(self, sigma=0.0, key=None):
        """ Returns all values where the chi-square is within the n-th sigma

        Parameters:
            sigma (float, int): Value of sigma to cut values.
            key (str): keyword the user desire to obtain results.

        Returns:
            - if a key is given, it returns list with all the values that are within the n-sigma.
            - if no key is given, it returns a dictionary with the list with all the values
                that are within the n-sigma for all keys.
            - if sigma is zero, it returns the parameters for the minimum chi-square instead of a list.
        """
        values = {}
        if sigma == 0.0:
            k = np.argsort(self.data['chi2'])[0]
        else:
            k = np.where(self.data['chi2'] < self.data['chi2'].min() + sigma**2)[0]
        for name in self.__names[1:]:
            values[name] = self.data[name][k]
        return values

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
        for name in self.__names[1:]:
            output += ('\n{}:\n'
                       '    1-sigma: {:.3f} +/- {:.3f}\n'
                       '    3-sigma: {:.3f} +/- {:.3f}\n'.format(
                           name, *sigma_1[name], *sigma_3[name])
                       )
        return output
