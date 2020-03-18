import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u

def draw_ellipse(radius, oblateness=0.0, center_x=0.0, center_y=0.0, pos_angle=0.0, *args, **kwargs):
    """ Plots an ellipse given input parameters

    Parameters:
        radius (float, int): Semi-major axis of the ellipse.
        oblateness (float, int): Oblateness of the ellipse. Default=0.0
        center_x (float, int): Coordinate of the ellipse (abscissa). Default=0.0
        center_y (float, int): Coordinate of the ellipse (ordinate). Default=0.0
        pos_angle (float, int): Pole position angle. Default=0.0
        **kwargs: all other parameters will be parsed directly to matplotlib
    """
    theta= np.linspace(-np.pi, np.pi, 1800)
    circle_x =  radius*np.cos(theta)
    circle_y =  radius*(1.0-oblateness)*np.sin(theta)
    if 'color' not in kwargs:
        kwargs['color'] = 'black'
    if 'lw' not in kwargs:
        kwargs['lw'] = 2
    plt.plot(circle_x*np.cos(pos_angle*u.deg) + circle_y*np.sin(pos_angle*u.deg) + center_x,
             -circle_x*np.sin(pos_angle*u.deg) + circle_y*np.cos(pos_angle*u.deg) + center_y,
             *args, **kwargs)
    
    
class ChiSquare():
    """ ChiSquare stores the arrays for all inputs and given chi-square.

    Parameters:
        chi2 (array): Array with all the chi-square values
        **kwargs: any other given input must be an array with the same size as chi2.
            the keyword name will be associated as the variable name of the given data
            
    Example:
    
    chisquare = ChiSquare(chi2, immersion=t1, emersion=t2)
    t1 and t2 must be an array with the same size as chi2.
    the data can be accessed as:
    chisquare.data['immersion']

    """
    def __init__(self,chi2,**kwargs):
        names = [('chi2', 'f8')]
        data = chi2
        data_size = len(chi2)
        for item in kwargs.keys():
            if len(kwargs[item]) != data_size:
                raise ValueError('{} size must have the same size as given chi2')
            names.append((item, 'f8'))
            data = np.vstack((data, kwargs[item]))
        self.data = np.array([tuple(i) for i in data.T], dtype=names)
        
    def get_nsigma(self,sigma=1, key=None):
        """ Determines the interval of the chi-square within the n-th sigma

        Parameters:
            sigma (float, int): Value of sigma to calculate.
            key (str): keyword the user desire to obtain results.
            
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
        for name in self.data.dtype.names[1:]:
            vmax = self.data[name][values].max()
            vmin = self.data[name][values].min()
            output[name] = [(vmax+vmin)/2.0, (vmax-vmin)/2.0]
        if key is not None:
            if key not in self.data.dtype.names[1:]:
                raise ValueError('{} is not one of the available keys. Please choose one of {}'
                                 .format(key, self.data.dtype.names[1:]))
            return output[key]
        return output
    
    def plot_chi2(self, key=None):
        """ Plots an ellipse given input parameters

        Parameters:
            key (str): Key to plot chi square.
            if no key  is given, it will plot for all the keywords.
        """
        sigma_1 = self.get_nsigma(sigma=1)
        sigma_3 = self.get_nsigma(sigma=3)
        if key is not None and (key not in self.data.dtype.names[1:]):
            raise ValueError('{} is not one of the available keys. Please choose one of {}'
                                 .format(key, self.data.dtype.names[1:]))
        for name in self.data.dtype.names[1:]:
            if (key is not None) and (key != name):
                continue
            plt.plot(self.data[name], self.data['chi2'], 'k.')
            plt.ylim(sigma_1['chi2_min']-1, sigma_1['chi2_min']+10)
            plt.xlim(sigma_3[name][0]-3*sigma_3[name][1], sigma_3[name][0]+3*sigma_3[name][1])
            plt.axhline(sigma_1['chi2_min']+1,linestyle='--',color='red')
            plt.axhline(sigma_3['chi2_min']+9,linestyle='--',color='red')
            #pl.axvline(t_i[chi2.argmin()],linestyle='--',color='red')
            plt.xlabel(name,fontsize=20)
            plt.ylabel('Chi2',fontsize=20)
            if key is None:
                plt.show()
                
    def to_file(self, namefile):
        """ Save the data to a file

        Parameters:
            namefile (str): Filename to save the data
        """
        np.savetxt(namefile, self.data, fmt='%10.5f')
        f = open(namefile+'.label', 'w')
        for i, name in enumerate(self.data.dtype.names):
            f.write('Column {}: {}\n'.format(i+1,name))
        f.close()
                
    def __str__(self):
        """ String representation of the ChiSquare Object
        """
        sigma_1 = self.get_nsigma(sigma=1)
        sigma_3 = self.get_nsigma(sigma=3)
        output = 'Minimum chi-square: {}\n'.format(sigma_1['chi2_min'])
        for name in self.data.dtype.names[1:]:
            output += '\n{}:\n'.format(name)
            output += '    1-sigma: {:.3f} +/- {:.3f}\n'.format(*sigma_1[name])
            output += '    3-sigma: {:.3f} +/- {:.3f}\n'.format(*sigma_3[name])
        return output