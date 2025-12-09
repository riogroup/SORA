import os
import warnings

import astropy.units as u
import numpy as np
from astropy.time import Time

from sora.config import input_tests
from sora.config.decorators import deprecated_alias
from .utils import calc_fresnel
from .occdetect import occ_detect
from .model import *
from .fit import fit
from .fit_double import fit_double

warnings.simplefilter('always', UserWarning)

class LightCurve:
    """ Defines a Light Curve.

    Parameters
    ----------
    name : `str`
        The name of the LightCurve. Each time a LightCurve object is defined
        the name must be different.

    tref : `astropy.time.Time`, `str`, `float`
        Instant of reference.

        Format: `Julian Date`, string in ISO format or Time object.
        Required only when input times are not given in Julian Date.

    central_bandpass : `int`, `float`, optional, default=0.58
        The effective central wavelength of the detector used in observation.
        Value in microns.

    delta_bandpass : `int`, `float`, optional, default=0.40
        The effective bandpass width of the detector used in observation.
        Value in microns.

    exptime : `int`, `float`
        The exposure time of the observation, in seconds.
        *Required* in case *1* below.
        *Not required* in cases *2*, *3* and *4* below.

    response : `tuple` of array_like, optional
        Instrumental or filter spectral response.
        Must be a tuple ``(lambda [µm], transmission)`` of the same length.

        When given, a Gauss-Legendre integration of the Fresnel diffraction is
        performed across the response bandpass. The number of nodes is set by
        `n_lambda` (minimum of 5). In this case, the uniform-bandpass
        approximation using ``delta_bandpass`` is ignored.

    n_lambda : `int`, optional, default=2
        Spectral sampling parameter controlling the number of Gauss-Legendre
        nodes used in the wavelength integration.

        If `response` is **None**:
            - ``n_lambda = 1`` → monochromatic model at `central_bandpass`.
            - ``n_lambda = 2`` → classic SORA two-endpoint approximation
              (``central_bandpass ± delta_bandpass/2``).
            - ``n_lambda ≥ 3`` → Gauss-Legendre integration over the top-hat
              bandpass.

        If `response` is **not None**:
            - Integration uses ``max(n_lambda, 5)`` nodes.
            - The physically accurate, response-weighted Fresnel diffraction
              is evaluated across the wavelength range defined by `response`.

    **kwargs : `int`, `float`
        Physical parameters of the occultation geometry.

        Note
        ----
        vel : `int`, `float`
            Shadow velocity, in km/s.

        dist : `int`, `float`
            Geocentric distance of the occulting body, in AU.

        d_star : `float`
            Stellar diameter, in km.


    Warning
    -------
    Input data must follow one of the 4 options below:

    1) Input data from file with time and flux  
        `file (str)`: Path to a file containing time and flux.  
        A third column with the flux error may also be given.

        `usecols (int, tuple or array)`: Which columns to read, with the first
        being time, the second flux, and the optional third the flux error.

        `skiprows (int)`: Number of header lines to skip.

        **Example:**

        >>> LightCurve(name, file=file_path, exptime=0.1)


    2) Input data when file is not given:  
        `time`: array-like of times (seconds from tref, JD, or Time objects).

        `flux`: array-like of fluxes. Must match the length of `time`.

        `dflux`: array-like of flux uncertainties (optional).

        **Example:**

        >>> LightCurve(name, time=time_array, flux=flux_array, exptime=0.1)


    Cases for when `time` and `flux` are not given.


    3) Input for a positive occultation:  
        `immersion`: The instant of immersion.  
        `emersion`: The instant of emersion.  
        `immersion_err`: Immersion time uncertainty, in seconds.  
        `emersion_err`: Emersion time uncertainty, in seconds.

        **Example:**

        >>> LightCurve(name, immersion, immersion_err, emersion, emersion_err)


    4) Input for a negative occultation:  
        `initial_time`: The initial time of observation.  
        `end_time`: The end time of observation.

        **Example:**

        >>> LightCurve(name, initial_time, end_time)


    Notes
    -----
    - If a `response` tuple is given, the effective bandpass (central_bandpass,
      delta_bandpass) is derived by integrating the normalized transmission curve.

    - Default filter parameters:
        central_bandpass = 0.58 µm  
        delta_bandpass   = 0.40 µm
    """   

    @deprecated_alias(lambda_0='central_bandpass', delta_lambda='delta_bandpass')  # remove this line for v1.0
    def __init__(self, name='', **kwargs):
        allowed_kwargs = ['immersion', 'immersion_err',
                          'emersion', 'emersion_err',  
                          'initial_time', 'end_time',
                          'file', 'time', 'flux', 
                          'exptime', 'dflux',
                          'central_bandpass', 'delta_bandpass', 
                          'response', 'n_lambda',
                          'tref', 'dist', 'vel', 'd_star', 
                          'skiprows', 'usecols']
        
        input_tests.check_kwargs(kwargs, allowed_kwargs=allowed_kwargs)
        input_done = False
        self.dflux = None
        self._name = name
        self.flux = None
        self.time_model = None
        if 'exptime' in kwargs:
            if kwargs['exptime'] <= 0:
                raise ValueError('Exposure time can not be zero or negative')
            self.exptime = kwargs['exptime']
        if 'vel' in kwargs:
            self.set_vel(vel=kwargs['vel'])
        if 'dist' in kwargs:
            self.set_dist(dist=kwargs['dist'])
        if 'd_star' in kwargs:
            self.set_star_diam(d_star=kwargs['d_star'])
        if 'tref' in kwargs:
            self.tref = kwargs['tref']
        if 'immersion' in kwargs:
            self.immersion = kwargs['immersion']
            self.immersion_err = kwargs.get('immersion_err', 0.0)
            if self.immersion_err < 0:
                warnings.warn("Immersion Error must be positive. Using absolute value.")
                self.immersion_err = np.absolute(self.immersion_err)
            input_done = True
        if 'emersion' in kwargs:
            self.emersion = kwargs['emersion']
            self.emersion_err = kwargs.get('emersion_err', 0.0)
            if self.emersion_err < 0:
                warnings.warn("Emersion Error must be positive. Using absolute value.")
                self.emersion_err = np.absolute(self.emersion_err)
            try:
                if self.emersion <= self.immersion:
                    raise ValueError("emersion time must be greater than immersion time")
            except AttributeError:
                pass
            input_done = True
        if 'initial_time' in kwargs and 'end_time' in kwargs:
            self.initial_time = kwargs['initial_time']
            self.end_time = kwargs['end_time']
            if self.end_time <= self.initial_time:
                raise ValueError('end_time must be greater than initial_time')
            input_done = True
        if 'flux' in kwargs or 'file' in kwargs:
            self.set_flux(**kwargs)
            input_done = True
        if not input_done:
            raise ValueError('No allowed input conditions satisfied. Please refer to the tutorial.')                
        if 'response' in kwargs:
            lam, T = kwargs['response']
            lam = np.asarray(lam, dtype=float)
            T = np.asarray(T, dtype=float)

            if len(lam) < 3 or np.all(T == 0):
                self.set_filter(
                    central_bandpass=kwargs.get('central_bandpass', 0.70),
                    delta_bandpass=kwargs.get('delta_bandpass', 0.30)
                )
            else:
                order = np.argsort(lam)
                lam, T = lam[order], T[order]
                norm = np.trapz(T, lam)
                if norm <= 0 or not np.isfinite(norm):
                    self.set_filter(
                        central_bandpass=kwargs.get('central_bandpass', 0.58),
                        delta_bandpass=kwargs.get('delta_bandpass', 0.40)
                    )
                else:
                    Tn = T / norm
                    lambda_0 = np.trapz(lam * Tn, lam)
                    sigma_lambda = np.sqrt(np.trapz((lam - lambda_0)**2 * Tn, lam))
                    delta_lambda = 2 * sigma_lambda
                    self.set_filter(central_bandpass=lambda_0, delta_bandpass=delta_lambda)
                self.set_response((lam, T))
        else:
            self.response = None
            self.set_filter(
                central_bandpass=kwargs.get('central_bandpass', 0.58),
                delta_bandpass=kwargs.get('delta_bandpass', 0.40)
            )

        self.dt = 0.0
        self.n_lambda = kwargs.get('n_lambda', 2)

    @property
    def fresnel_scale(self):
        lamb = self.lambda_0*u.micrometer.to('km')
        dlamb = self.delta_lambda*u.micrometer.to('km')
        dist = self.dist*u.au.to('km')
        fresnel_scale_1 = calc_fresnel(dist, lamb-dlamb/2.0)
        fresnel_scale_2 = calc_fresnel(dist, lamb+dlamb/2.0)
        fresnel_scale = (fresnel_scale_1 + fresnel_scale_2)/2.0
        return fresnel_scale

    @property
    def central_bandpass(self):
        return self.lambda_0

    @property
    def delta_bandpass(self):
        return self.delta_lambda

    @property
    def name(self):
        return self._name

    @property
    def tref(self):
        if hasattr(self, '_tref'):
            return self._tref
        else:
            raise AttributeError("'LightCurve' object has no attribute 'tref'")

    @tref.setter
    def tref(self, value):
        if type(value) in [int, float]:
            self.tref = Time(value, format='jd')
        else:
            try:
                self._tref = Time(value)
            except ValueError:
                raise ValueError('{} is not a valid time format accepted by tref'.format(value))

    @property
    def immersion(self):
        if hasattr(self, '_immersion'):
            return self._immersion + self.dt*u.s
        else:
            raise AttributeError('The immersion time was not fitted or instantiated.')

    @immersion.setter
    def immersion(self, value):
        if type(value) in [int, float]:
            if value > 2400000:
                self._immersion = Time(value, format='jd')
            elif hasattr(self, 'tref'):
                self._immersion = self.tref + value*u.s
            else:
                raise ValueError('{} can not be set without a reference time'.format(value))
        else:
            try:
                self._immersion = Time(value)
            except ValueError:
                raise ValueError('{} is not a valid time format accepted by immersion'.format(value))

    @property
    def emersion(self):
        if hasattr(self, '_emersion'):
            return self._emersion + self.dt*u.s
        else:
            raise AttributeError('The emersion time was not fitted or instanciated.')

    @emersion.setter
    def emersion(self, value):
        if type(value) in [int, float]:
            if value > 2400000:
                self._emersion = Time(value, format='jd')
            elif hasattr(self, 'tref'):
                self._emersion = self.tref + value*u.s
            else:
                raise ValueError('{} can not be set without a reference time'.format(value))
        else:
            try:
                self._emersion = Time(value)
            except ValueError:
                raise ValueError('{} is not a valid time format accepted by emersion'.format(value))

    @property
    def initial_time(self):
        if hasattr(self, '_initial_time'):
            return self._initial_time
        else:
            raise AttributeError("'LightCurve' object has no attribute 'initial_time'")

    @initial_time.setter
    def initial_time(self, value):
        if type(value) in [int, float]:
            if value > 2400000:
                self._initial_time = Time(value, format='jd')
            elif hasattr(self, 'tref'):
                self._initial_time = self.tref + value*u.s
            else:
                raise ValueError('{} can not be set without a reference time'.format(value))
        else:
            try:
                self._initial_time = Time(value)
            except ValueError:
                raise ValueError('{} is not a valid time format accepted by initial_time'.format(value))

    @property
    def end_time(self):
        if hasattr(self, '_end_time'):
            return self._end_time
        else:
            raise AttributeError("'LightCurve' object has no attribute 'end_time'")

    @end_time.setter
    def end_time(self, value):
        if type(value) in [int, float]:
            if value > 2400000:
                self._end_time = Time(value, format='jd')
            elif hasattr(self, 'tref'):
                self._end_time = self.tref + value*u.s
            else:
                raise ValueError('{} can not be set without a reference time'.format(value))
        else:
            try:
                self._end_time = Time(value)
            except ValueError:
                raise ValueError('{} is not a valid time format accepted by end_time'.format(value))

    @property
    def time_mean(self):
        if hasattr(self, '_immersion') and hasattr(self, '_emersion'):
            return Time((self.immersion.jd + self.emersion.jd)/2, format='jd')
        else:
            return Time((self.initial_time.jd + self.end_time.jd)/2, format='jd')

    @property
    def time(self):
        try:
            return (self._time - self.tref).sec + self.dt
        except:
            raise AttributeError("'LightCurve' object has no attribute 'time'")
        
    @property
    def detection_limits(self):
        if not hasattr(self, "_detection_limits"):
            from sora.rings.detection_limits import DetectionLimits
            self._detection_limits = DetectionLimits(self)
        return self._detection_limits

    def set_flux(self, **kwargs):
        """ Sets the flux for the LightCurve.

        Parameters
        ----------
        exptime : `int`, `float`, required
            The exposure time of the observation, in seconds.

        file : `str`
            A file with the time and flux in the first and second columns,
            respectively. A third column with error in flux can also be given.

        time
            If file not given, time must be a list of times, in seconds from `tref`,
            or `Julian Date`, or a `Time object`.

        flux
            If file not given, flux must be a list of fluxes. It must have the
            same lenght as time.

        dflux
            If file not given, dflux must be a list of fluxes errors. It must
            have the same lenght as time.

        tref : `astropy.time.Time`, `str`, `float`
            Instant of reference. It can be in `Julian Date`, string in ISO
            format or `Time object`.

        usecols : `int`, `tuple`, array, optional
            Which columns to read, with the first being the time, the seconds
            the flux and third the flux error.

        **kwargs : `int`, `float`
            Object velocity, object distance, star diameter.

            Note
            ----
            vel : `int`, `float`
                Velocity in km/s.

            dist : `int`, `float`:
                Object distance in AU.

            d_star : `float`
                Star diameter, in km.
        """
        from .utils import read_lc_file
        import scipy.stats as scst

        input_done = False
        usecols = None
        if 'usecols' in kwargs:
            usecols = kwargs['usecols']
        skiprows = 0
        if 'skiprows' in kwargs:
            skiprows = int(kwargs['skiprows'])
        if 'file' in kwargs:
            try:
                lc_data = read_lc_file(kwargs['file'], usecols=usecols, skiprows=skiprows)
                if len(lc_data) == 2:
                    time, self.flux = lc_data
                elif len(lc_data) == 3:
                    time, self.flux, self.dflux = lc_data
            except:
                pass
            if hasattr(self, 'flux'):
                self.flux_obs = self.flux
            if not hasattr(self, 'flux_obs'):
                raise ValueError('Input file must have 2 or 3 columns')
            input_done = True
        if 'time' in kwargs and 'flux' in kwargs:
            if input_done:
                raise ValueError('Only one type of input can be given. Please refer to the tutorial.')
            self.flux = kwargs['flux']
            time = kwargs['time']
            if len(self.flux) != len(time):
                raise ValueError('time and flux must have the same length')
            if 'dflux' in kwargs:
                self.dflux = kwargs['dflux']
                if len(self.flux) != len(self.dflux):
                    raise ValueError('dflux must have the same length as flux and time')
            input_done = True
        if not input_done:
            raise ValueError('Input parameters not satisfied')
        if 'exptime' not in kwargs:
            raise ValueError('exptime not defined')
        if kwargs['exptime'] <= 0:
            raise ValueError('Exposure time can not be zero or negative')
        else:
            self.exptime = kwargs['exptime']
        if 'vel' in kwargs:
            self.set_vel(vel=kwargs['vel'])
        if 'dist' in kwargs:
            self.set_dist(dist=kwargs['dist'])
        if 'd_star' in kwargs:
            self.set_star_diam(d_star=kwargs['d_star'])
        if 'tref' in kwargs:
            self.tref = kwargs['tref']
        if 'time' in locals():
            if type(time) == Time:
                if not hasattr(self, 'tref'):
                    self.tref = Time(time[0].iso.split(' ')[0] + ' 00:00:00.000')
            elif all(time > 2400000):
                time = Time(time, format='jd')
                if not hasattr(self, 'tref'):
                    self.tref = Time(time[0].iso.split(' ')[0] + ' 00:00:00.000')
            elif not hasattr(self, 'tref'):
                raise ValueError('tref must be given')
            else:
                time = self.tref + time*u.s
            order = np.argsort(time)
            self._time = time[order]
            self.flux = self.flux[order]
            self.flux_obs = self.flux
            if self.dflux is not None:
                self.dflux = self.dflux[order]
            self.initial_time = np.min(time)
            self.end_time = np.max(time)
            time_diffs = time_diffs = (self._time[1:] - self._time[:-1]).sec
            self.cycle = scst.mode(time_diffs, keepdims=False).mode
            if self.cycle < self.exptime:
                warnings.warn('Exposure time ({:0.4f} seconds) higher than Cycle time ({:0.4f} seconds)'.
                              format(self.exptime, self.cycle))

    def set_exptime(self, exptime):
        """ Sets the light curve exposure time.

        Parameters
        ----------
        exptime : `int`, `float`
            Exposure time, in seconds.
        """
        exptime = u.Quantity(exptime, unit=u.s)
        if not np.isscalar(exptime):
            raise TypeError('Exposure time must be an integer, a float or an Astropy Unit object')
        if exptime.value <= 0:
            raise ValueError('Exposure time can not be zero or negative')
        self.exptime = exptime.value
        try:
            if self.cycle < self.exptime:
                warnings.warn('Exposure time ({:0.4f} seconds) higher than Cycle time ({:0.4f} seconds)'.
                              format(self.exptime, self.cycle))
        except:
            pass

    def set_vel(self, vel):
        """ Sets the occultation velocity.

        Parameters
        ----------
        vel : `int`, `float`
            Velocity in km/s.
        """
        vel = u.Quantity(vel, unit=u.km/u.s)
        self.vel = np.absolute(vel.value)

    def set_dist(self, dist):
        """ Sets the object distance.

        Parameters
        ----------
        dist : `int`, `float`
            Object distance in AU.
        """
        dist = u.Quantity(dist, unit=u.AU)
        if dist.value < 0:
            warnings.warn("distance cannot be negative. Using absolute value.")
        self.dist = np.absolute(dist.value)

    def set_star_diam(self, d_star):
        """ Sets the star diameter.

        Parameters
        ----------
        d_star : `float`
            Star diameter, in km.
        """
        d_star = u.Quantity(d_star, unit=u.km)
        if d_star.value < 0:
            warnings.warn("star diameter cannot be negative. Using absolute value.")
        self.d_star = np.absolute(d_star.value)

    @deprecated_alias(lambda_0='central_bandpass', delta_lambda='delta_bandpass')  # remove this line for v1.0
    def set_filter(self, central_bandpass, delta_bandpass):
        """ Sets the filter bandwidth in microns.

        Parameters
        ----------
        central_bandpass : `float`
            Center band in microns.

        delta_bandpass : `float`
            Bandwidth in microns.
        """
        central_bandpass = u.Quantity(central_bandpass, unit=u.micrometer)
        if central_bandpass.value <= 0:
            raise ValueError("central bandpass cannot be negative.")
        self.lambda_0 = central_bandpass.value
        delta_bandpass = u.Quantity(delta_bandpass, unit=u.micrometer)
        if delta_bandpass <= 0:
            raise ValueError("delta bandpass cannot be negative")
        self.delta_lambda = delta_bandpass.value
        if (central_bandpass - delta_bandpass).value <= 0:
            raise ValueError("The given central and delta bandpass give a range ({}, {}) microns. Bandpass cannot be negative. "
                             "Please give appropriate values".format(*(central_bandpass +
                                                                       np.array([-1, 1])*delta_bandpass).value))

    def set_response(self, response):
        """ Sets the spectral response curve (λ, Tλ).

        Parameters
        ----------
        response : array-like, tuple
            Two arrays (λ, Tλ) in microns and relative transmission (0-1).
        """
        import numpy as np

        lam, trans = response
        lam = np.asarray(lam, dtype=float)
        trans = np.asarray(trans, dtype=float)

        if lam.ndim != 1 or trans.ndim != 1 or lam.size != trans.size:
            raise ValueError("Response must be two 1D arrays of equal length (λ, Tλ).")

        if np.any(lam <= 0):
            raise ValueError("Wavelengths must be positive (μm).")

        if np.any(trans < 0) or np.any(trans > 1.2):
            raise ValueError("Transmission values must be between 0 and 1.")

        trans /= np.nanmax(trans)
        self.response = (lam, trans)

    def calc_magnitude_drop(self, mag_star, mag_obj):
        """ Determines the magnitude drop of the occultation.

        Parameters
        ----------
        mag_star : `int`, `float`
            Star magnitude.

        mag_obj `int`, `float`
            Object apparent magnitude to the date.

        Returns
        -------
        mag_drop : `float`
            Magnitude drop for the given magnitudes.

        bottom_flux : `float`
            Normalized bottom flux for the given magnitudes.
        """
        from .utils import calc_magnitude_drop
        mag_drop, bottom_flux = calc_magnitude_drop(mag_star, mag_obj)
        self.mag_drop = mag_drop
        self.bottom_flux = bottom_flux

    def normalize(self, poly_deg=None, mask=None, flux_min=0.0, flux_max=1.0, plot=False):
        """ Returns the normalized flux within the flux min and flux max defined scale.

        Parameters
        ----------
        poly_deg : `int`
            Degree of the polynomial to be fitted.

        mask : `bool` array
            Which values to be fitted.

        flux_min : `int`, `float`
            Event flux to be set as 0.

        flux_max : `int`, `float`
            Baseline flux to be set as 1.

        plot : `bool`
            If True plot the steps for visual aid.
        """
        from .utils import fit_pol
        import matplotlib.pyplot as plt

        if type(self.flux) == type(None):
            raise ValueError('Normalization is only possible when a LightCurve is instantiated with time and flux.')
        self.reset_flux()
        lc_flux = (self.flux - flux_min)/(flux_max-flux_min)
        if mask is None:
            preliminar_occ = occ_detect(self.flux, self.dflux, self.time, self.cycle,
                                    maximum_duration=((self.end_time - self.initial_time).to(u.s).value)/3)        
            tmax = preliminar_occ['emersion_time']+1.00*preliminar_occ['occultation_duration']
            tmin = preliminar_occ['immersion_time']-1.00*preliminar_occ['occultation_duration']
            chord = preliminar_occ['occultation_duration']
            mask = np.invert((self.time > tmin-(chord/2)) & (self.time < tmax+(chord/2)))
        norm_time = (self.time - self.time.min())/(self.time.max()-self.time.min())
        if poly_deg is not None:
            n = poly_deg
            p, err = fit_pol(norm_time[mask], lc_flux[mask], n)
            flux_poly_model = np.zeros(len(norm_time))
            for ii in np.arange(n+1):
                flux_poly_model = flux_poly_model + p[ii]*(norm_time**(n-ii))
            if plot:
                plt.plot(norm_time[mask], lc_flux[mask], 'k.-')
                plt.plot(norm_time[mask], flux_poly_model[mask], 'r-')
                plt.title('Polynomial degree = {}'.format(n), fontsize=15)
                plt.show()
        else:
            n = 0
            p, err = fit_pol(norm_time[mask], lc_flux[mask], n)
            flux_poly_model = np.zeros(len(norm_time))
            for ii in np.arange(n+1):
                flux_poly_model += p[ii]*(norm_time**(n-ii))
            if plot:
                plt.plot(norm_time[mask], lc_flux[mask], 'k.-')
                plt.plot(norm_time[mask], flux_poly_model[mask], 'r-')
                plt.title('Polynomial degree = {}'.format(n), fontsize=15)
                plt.show()
            for nn in np.arange(1, 10):
                p, err = fit_pol(norm_time[mask], lc_flux[mask], nn)
                flux_poly_model_new = np.zeros(len(norm_time))
                for ii in np.arange(nn+1):
                    flux_poly_model_new += p[ii]*(norm_time**(nn-ii))
                F = np.var(flux_poly_model[mask]-lc_flux[mask])/np.var(flux_poly_model_new[mask]-lc_flux[mask])
                if F > 1.05:
                    flux_poly_model = flux_poly_model_new.copy()
                    n = nn
                    if plot:
                        plt.plot(norm_time[mask], lc_flux[mask], 'k.-')
                        plt.plot(norm_time[mask], flux_poly_model[mask], 'r-')
                        plt.title('Polynomial degree = {}'.format(nn), fontsize=15)
                        plt.show()
                else:
                    print('Normalization using a {} degree polynomial'.format(n))
                    print('There is no improvement with a {} degree polynomial'.format(n+1))
                    break
        self.flux = lc_flux/flux_poly_model
        self.normalizer_flux = flux_poly_model
        self.normalizer_mask = mask

    def reset_flux(self):
        """ Resets flux for original values
        """
        try:
            self.flux = self.flux_obs
        except:
            raise ValueError('Reset is only possible when a LightCurve is instantiated with time and flux.')
        return

    def occ_detect(self, maximum_duration=None, dur_step=None, snr_limit=None,
               n_detections=None, tmin=None, tmax=None, plot=False):
        """ Automatically detect an occultation event in the light curve.

        This is used internally by `fit()` to estimate immersion/emersion times
        when not provided by the user.

        Parameters
        ----------
        maximum_duration : float, optional
            Maximum duration (s) to consider for occultation.
        dur_step : float, optional
            Step size (s) for duration scan.
        snr_limit : float, optional
            Minimum signal-to-noise ratio.
        n_detections : int, optional
            Maximum number of detections to return.
        tmin, tmax : float, optional
            Time limits for search.
        plot : bool, optional
            If True, shows diagnostic plot.

        Returns
        -------
        dict
            Ordered dictionary with parameters:
            'immersion_time', 'emersion_time', 'depth', 'baseline',
            'occultation_duration', 'time_err', 'snr', etc.
        """
        from .occdetect import occ_detect

        occ = occ_detect(
            self.flux,
            self.dflux,
            self.time,
            getattr(self, "cycle", None),
            maximum_duration=maximum_duration,
            dur_step=dur_step,
            snr_limit=snr_limit,
            n_detections=n_detections,
            tmin=tmin,
            tmax=tmax,
            plot=plot,
        )
        return occ

    def plot(self, ax=None):
        """ Plots the light curve
        """
        import matplotlib.pyplot as plt

        if not hasattr(self, 'flux') or self.flux is None:
            raise ValueError("LightCurve must have time and flux defined to plot.")

        ax = ax or plt.gca()

        ax.plot(self.time, self.flux, 'k.-', label='Observed', zorder=0)
        
        if not hasattr(self, 'model') or not self.model.any():
            ax.plot(self.time, np.ones(len(self.time)), 'r-', label='Model', zorder=0)
            ax.scatter(self.time, np.ones(len(self.time)), edgecolor='r', facecolor='none', zorder=0)   
        else:
            ax.plot(self.time, self.model, 'r.-', label='Model', zorder=0)
            ax.scatter(self.time, self.model, edgecolor='r', facecolor='none', zorder=0) 
            
        ax.set_xlabel('Time [seconds]', fontsize=20)
        ax.set_ylabel('Relative Flux', fontsize=20)
        ax.set_title(self.name)
        ax.legend()


    def occ_model(self, immersion, emersion, opacity,
                npt_star=12, time_resolution_factor=10, flux_min=0, flux_max=1):
        """ Deprecated method (use `lc.SquareWellModel()` instead).

        Returns the modelled light curve using the Fresnel-aware square-well model.

        Parameters
        ----------
        immersion, emersion : float
            Immersion and emersion times (seconds relative to tref).
        opacity : float
            Opacity of the occulting region (1=opaque, 0=transparent).
        mask : array(bool)
            Boolean mask of the light curve region to model.
        npt_star, time_resolution_factor, flux_min, flux_max : optional
            Model parameters identical to the old implementation.

        Returns
        -------
        flux : ndarray
            Modeled flux array matching the input mask region.
        """
        warnings.warn(
            "LightCurve.occ_model() is deprecated and will be removed in a future release.\n"
            "Use `lc.SquareWellModel(immersion_time=..., emersion_time=..., opacity=...).compute(lc.time)` instead.",
            DeprecationWarning,
            stacklevel=2
        )

        model = SquareWellModel(
            lightcurve=self,
            immersion=immersion,
            emersion=emersion,
            opacity=opacity,
            npt_star=npt_star,
            time_resolution_factor=time_resolution_factor,
            flux_min=flux_min,
            flux_max=flux_max,
        )
        self.model = model.compute()
        return model
    
    def clear_fits(self):
        """ Remove all stored models / chi2 / fit results."""
        if hasattr(self, "models"):
            self.models.clear()
        if hasattr(self, "chi2_maps"):
            self.chi2_maps.clear()
        if hasattr(self, "_fit_results"):
            self._fit_results.clear()

    def remove_fit(self, label: str):
        """ Remove a specific fit labelled 'fit1', 'fit2', etc.
        Check the labelled fits using lightcurve.models
        """
        for attr in ("models", "chi2_maps", "_fit_results"):
            d = getattr(self, attr, None)
            if isinstance(d, dict) and label in d:
                del d[label]

    def to_log(self, namefile=None):
        """ Saves the light curve log to a file.

        Parameters
        ----------
        namefile : `str`
            Filename to save the log.
        """
        if namefile is None:
            namefile = self.name.replace(' ', '_')+'.log'
        f = open(namefile, 'w')
        f.write(self.__str__())
        f.close()

    def to_file(self, namefile=None, overwrite=False):
        """ Saves the observed light curve to an ASCII file with metadata.
        
        """
        if self.flux is None:
            raise ValueError("Cannot save light curve — no flux data found.")

        if namefile is None:
            namefile = f'{self.tref.iso[:10].replace("-", "")}_{self.name.replace(" ", "_").replace("-", "_")}_lightcurve.dat'

        if os.path.exists(namefile) and not overwrite:
            raise FileExistsError(f"File '{namefile}' already exists. Use overwrite=True to replace it.")

        header_lines = [f"SORA LightCurve export",
                        f"Light curve name: {self.name}"]

        if hasattr(self, "tref"):
            header_lines.append(f"Reference time (UTC): {self.tref.iso}")
        if hasattr(self, "initial_time") and hasattr(self, "end_time"):
            duration = (self.end_time - self.initial_time).to(u.min).value
            header_lines.append(f"Observation start: {self.initial_time.iso}")
            header_lines.append(f"Observation end:   {self.end_time.iso}")
            header_lines.append(f"Duration:          {duration:.2f} minutes")

        if hasattr(self, "exptime"):
            header_lines.append(f"Exposure time:     {self.exptime:.3f} s")
        if hasattr(self, "cycle"):
            header_lines.append(f"Cycle time:        {self.cycle:.3f} s")

        if hasattr(self, "vel"):
            header_lines.append(f"Shadow velocity:   {self.vel:.3f} km/s")
        if hasattr(self, "dist"):
            header_lines.append(f"Object distance:   {self.dist:.3f} AU")

        if hasattr(self, "lambda_0") and hasattr(self, "delta_lambda"):
            header_lines.append(f"Bandpass:          {self.lambda_0:.2f} ± {self.delta_lambda:.2f} μm")

        if hasattr(self, "d_star"):
            header_lines.append(f"Stellar diameter:  {self.d_star:.3f} km")

        header_lines.append("")
        header_lines.append("Columns: jd, sec from tref, flux, flux_uncertainty, modeled flux, residuals")
        header = "\n".join(header_lines)

        time_sec = self.time
        time_iso = Time(self.tref) + self.time*u.s
        time_jd = time_iso.jd
        for model in self.models:
            mdl_flux = self.models[model].model_flux    
            
        data = np.column_stack((time_jd, 
                                time_sec, 
                                self.flux, 
                                self.dflux if self.dflux is not None else np.repeat(self.flux.std(ddof=1), len(self.flux)),
                                mdl_flux, self.flux - mdl_flux))

        np.savetxt(namefile, data, header=header, fmt="%.6f")

        return namefile

    def __str__(self):
        """ String representation of the LightCurve Object
        """
        output = 'Light curve name: {}\n'.format(self.name)
        try:
            output += ('Initial time: {} UTC\n'
                    'End time:     {} UTC\n'
                    'Duration:     {:.3f} minutes\n'.format(
                        self.initial_time.iso, self.end_time.iso,
                        (self.end_time - self.initial_time).value*u.d.to('min'))
                    )
        except:
            pass
        output += 'Time offset:  {:.3f} seconds\n\n'.format(self.dt)
        try:
            output += 'Exposure time:    {:.4f} seconds\n'.format(self.exptime)
            output += 'Cycle time:       {:.4f} seconds\n'.format(self.cycle)
            output += 'Num. data points: {}\n\n'.format(len(self.time))
        except:
            output += 'Object LightCurve was not instantiated with time and flux.\n\n'

        try:
            if getattr(self, "response", None) is not None:
                lam, T = self.response
                output += (
                    f"Spectral response:    Custom curve ({len(lam)} points)\n"
                    f"  Bandpass (derived): {self.lambda_0:.3f} ± {self.delta_lambda:.3f} μm\n"
                )
                if self.n_lambda != 2:   
                    output += f"  λ-sampling:         Gauss-Legendre, N = n_lambda = {self.n_lambda}\n"
            else:
                output += (
                    f"Bandpass:             {self.lambda_0:.3f} ± {self.delta_lambda:.3f} μm\n"
                )
                if self.n_lambda != 2:   
                    output += f"  λ-sampling:         n_lambda = {self.n_lambda}\n"


            output += (
                f"Object Distance:      {self.dist:.2f} AU\n"
                f"Used shadow velocity: {self.vel:.3f} km/s\n"
                f"Fresnel scale:        {self.fresnel_scale/self.vel:.3f} seconds or {self.fresnel_scale:.2f} km\n"
                f"Stellar size effect:  {self.d_star/self.vel:.3f} seconds or {self.d_star:.2f} km\n"
            )
        except Exception as e:
            output += f'\nThere is no occultation associated with this light curve {e}.\n'

        try:
            lc_dl = self.detection_limits
            output += lc_dl.__str__()
        except:
            pass
                    
        if hasattr(self, "_fit_results") and len(self._fit_results) > 0:
            output += "\nFitted models:\n"
            for label, res in self._fit_results.items():
                kind = res.get("type", "SquareWell")
                output += f"  {label} ({kind}):\n"
                try:
                    output += ('    Inst. response:       {:.3f} seconds or {:.2f} km\n'
                               '    Dead time effect:     {:.3f} seconds or {:.2f} km\n'
                               #'    Model resolution:     {:.3f} seconds or {:.2f} km\n'
                               '    Modelled baseflux:    {:.3f}\n'
                               '    Modelled bottomflux:  {:.3f}\n'
                               '    Light curve sigma:    {:.3f}\n\n'.format(
                                self.exptime, self.exptime*self.vel, 
                                self.cycle-self.exptime, (self.cycle-self.exptime)*self.vel, 
                                #self.model_resolution, self.model_resolution*self.vel, 
                                res['baseflux'], 
                                res['bottomflux'],
                                res['curve_sigma'])
                            )
                except Exception as e:
                    print(f"[DEBUG] LightCurve __str__ format error: {e}")

                if kind == "DoubleSquareWell":
                    output += (f"    Immersion 1: {res['immersion1_time'].iso}  Op1={res['opacity1']:.3f}\n"
                            f"    Emersion 1:  {res['emersion1_time'].iso}\n"
                            f"    Immersion 2: {res['immersion2_time'].iso}  Op2={res['opacity2']:.3f}\n"
                            f"    Emersion 2:  {res['emersion2_time'].iso}\n\n")
                else:
                    output += (f"    Immersion: {res['immersion_time'].iso} ± {res.get('immersion_err', 0):.3f} s\n"
                            f"    Emersion:  {res['emersion_time'].iso} ± {res.get('emersion_err', 0):.3f} s\n"
                            f"    Opacity:   {res['opacity']:.3f}\n\n")       
                try:
                    if hasattr(self, "chi2_maps") and label in self.chi2_maps:
                        chi2 = self.chi2_maps[label]
                        try:
                            chi_text = chi2.__str__().rstrip()
                            chi_text = "\n".join("    " + line for line in chi_text.splitlines())
                            output += chi_text + "\n\n"
                        except Exception:
                            pass
                except Exception:
                    pass
        else:
            output += '\nObject LightCurve model was not fitted.\n\n'

        return output 

LightCurve.fit = fit
LightCurve.fit_double = fit_double
attach_to_lightcurve_class(LightCurve)