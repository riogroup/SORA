import os

import numpy as np

from sora.config.decorators import deprecated_alias

__all__ = ['calc_fresnel', 'calc_magnitude_drop']


@deprecated_alias(lambida='bandpass')  # remove this line for v1.0
def calc_fresnel(distance, bandpass):
    """Calculates the Fresnel scale.

    Fresnel Scale = square root of half the multiplication of wavelength and
    object distance.

    Parameters
    ----------
    distance : `int`, `float` array
        Distances, in km.

    bandpass : `int`, `float`, array
        Wavelength, in km.

    Returns
    -------
    fresnel_scale : `float`, array
        Fresnel scale, in km.
    """
    return np.sqrt(bandpass * distance / 2)


def bar_fresnel(X, X01, X02, fresnel_scale, opacity):
    """Returns the modelled light curve considering fresnel diffraction.

    Parameters
    ----------
    X : array
        Array with time values converted in km using the event velocity.

    X01 : `int`, `float`
        Immersion time converted in km using the event velocity.

    X02 `int`, `float`
        Emersion time converted in km using the event velocity.

    fresnel_scale : `int`, `float`
        Fresnel scale, in km.

    opacity : `int`, `float`
        Opacity. Opaque = 1.0, transparent = 0.0.

    Returns
    -------
    flux_fresnel : array
        The light curve with fresnel diffraction.
    """
    import scipy.special as scsp

    # Converting from km to units of fresnel scale
    x = X / fresnel_scale
    x01 = X01 / fresnel_scale
    x02 = X02 / fresnel_scale
    # Fresnel diffraction parameters
    x1 = x - x01
    x2 = x - x02
    s1, c1 = scsp.fresnel(x1)
    s2, c2 = scsp.fresnel(x2)
    cc = c1 - c2
    ss = s1 - s2
    r_ampli = - (cc + ss) * (opacity / 2.)
    i_ampli = (cc - ss) * (opacity / 2.)
    # Determining the flux considering fresnel diffraction
    flux_fresnel = (1.0 + r_ampli) ** 2 + i_ampli ** 2
    return flux_fresnel


def fit_pol(x, y, deg):
    """Fits a n-degree polynomial to the data.

    Parameters
    ----------
    x : `array`
        x-values.

    y : `array`
        y-values.

    deg : `int`
        Degree of the polynomial.


    Returns
    -------
    param : array
        Array with the fitted values.

    param_err : array
        Array with the errors of the fitted values.
    """
    from scipy.odr import odrpack as odr
    from scipy.odr import models

    func = models.polynomial(deg)
    mydata = odr.Data(x, y)
    myodr = odr.ODR(mydata, func, maxit=200)
    myodr.set_job(fit_type=2)
    fit = myodr.run()
    param = fit.beta[::-1]
    param_err = fit.sd_beta[::-1]
    return param, param_err


def calc_magnitude_drop(mag_star, mag_obj):
    """Determines the magnitude drop of the occultation.

    Parameters
    ----------
    mag_star : `int`, `float`
        Star magnitude.

    mag_obj : `int`, `float`
        Object apparent magnitude to the date.


    Returns
    -------
    mag_drop : `float`
        Magnitude drop for the given magnitudes.

    bottom_flux : `float`
        Normalized bottom flux for the given magnitudes.
    """
    contrast = 1 / (1 + (10 ** ((mag_star - mag_obj) * 0.4)))
    mag_combined = mag_star - (2.5 * (np.log10(1 / contrast)))
    mag_drop = mag_obj - mag_combined
    bottom_flux = 10 ** ((mag_combined - mag_obj) * 0.4)
    mag_drop = mag_drop
    bottom_flux = bottom_flux
    return mag_drop, bottom_flux


def read_lc_file(lc_file, usecols=None, skiprows=0):
    """Reads light curve data from a file.

    Parameters
    ----------
    file : `str`, required
        A file with the time and flux in the first and second columns, respectively.
        A third column with error in flux can also be given.

    usecols : `int`, `tuple`, array, optional
        Which columns to read, with the first being the time, the seconds the
        flux and third the flux error.

    skiprows : `int`, optional, default=0
        Skip the first skiprows lines, including comments.

    Returns
    -------
    :rtype: : (array, array, array)

    """
    from astropy.io import ascii
    from astropy.time import Time

    if not os.path.isfile(lc_file):
        raise ValueError('{} not found'.format(lc_file))
    if usecols is not None:
        if len(usecols) == 2:
            try:
                time, flux = np.loadtxt(lc_file, usecols=usecols, skiprows=skiprows, unpack=True)
                return time, flux
            except:
                pass
        elif len(usecols) == 3:
            try:
                time, flux, dflux = np.loadtxt(lc_file, usecols=usecols, skiprows=skiprows, unpack=True)
                return time, flux, dflux
            except:
                pass
        else:
            raise ValueError('usecols should have 2 or 3 values')
    else:
        try:
            time, flux = np.loadtxt(lc_file, usecols=[0, 1], skiprows=skiprows, unpack=True)
            return time, flux
        except:
            pass
        try:
            time, flux, dflux = np.loadtxt(lc_file, usecols=[0, 1, 2], skiprows=skiprows, unpack=True)
            return time, flux, dflux
        except:
            pass
    try:
        lc_data = ascii.read(lc_file, data_start=skiprows)

        if usecols is not None:
            if len(usecols) == 2:
                time_col_index, flux_col_index = usecols
                time_col_name = lc_data.colnames[time_col_index]
                flux_col_name = lc_data.colnames[flux_col_index]
                time, flux = Time(lc_data[time_col_name].data).jd, lc_data[flux_col_name].data
                return time, flux
            elif len(usecols) == 3:
                time_col_index, flux_col_index, dflux_col_index = usecols
                time_col_name = lc_data.colnames[time_col_index]
                flux_col_name = lc_data.colnames[flux_col_index]
                dflux_col_name = lc_data.colnames[dflux_col_index]
                time, flux = Time(lc_data[time_col_name].data).jd, lc_data[flux_col_name].data
                dflux = lc_data[dflux_col_name].data
                return time, flux, dflux
            else:
                raise ValueError('usecols should have 2 or 3 values')
        else:
            if len(lc_data.colnames) == 2:
                time_col_index, flux_col_index = [0, 1]
                time_col_name = lc_data.colnames[time_col_index]
                flux_col_name = lc_data.colnames[flux_col_index]
                time, flux = Time(lc_data[time_col_name].data).jd, lc_data[flux_col_name].data
                return time, flux
            elif len(lc_data.colnames) == 3:
                time_col_index, flux_col_index, dflux_col_index = [0, 1, 2]
                time_col_name = lc_data.colnames[time_col_index]
                flux_col_name = lc_data.colnames[flux_col_index]
                dflux_col_name = lc_data.colnames[dflux_col_index]
                time, flux = Time(lc_data[time_col_name].data).jd, lc_data[flux_col_name].data
                dflux = lc_data[dflux_col_name].data
                return time, flux, dflux
            else:
                raise ValueError('file columns  should have 2 or 3 values')
    except:
        raise ValueError('An error has occurred in reading the {}'.format(lc_file))
