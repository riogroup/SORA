import os
import warnings

import astropy.units as u
import numpy as np
from astropy.time import Time

from sora.config import input_tests
from sora.config.decorators import deprecated_alias
from .utils import calc_fresnel

warnings.simplefilter('always', UserWarning)


class LightCurve:
    """Defines a Light Curve.

    Parameters
    ----------
    name : `str`
        The name of the LightCurve. Each time an LightCurve object is defined
        the name must be different.

    tref : `astropy.time.Time`, `str`, `float`
        Instant of reference.

        Format: `Julian Date`, string in ISO format or Time object.
        Required only if LightCurve have input fluxes and given time is
        not in Julian Date.

    central_bandpass : `int`, `float`, otpional, default=0.7
        The center band pass of the detector used in observation. Value in microns.

    delta_bandpass : `int`, `float`, optional, default=0.3
        The band pass width of the detector used in observation. Value in microns.

    exptime : `int`, `float`
        The exposure time of the observation, in seconds.
        *NOT* required in cases *2*, *3* and *4* below.
        *Required* in case *1* below.

    **kwargs: `int`, `float`
        Object velocity, distance, and star diameter.

        Note
        ----
        vel : `int`, `float`
            Velocity in km/s.

        dist : `int`, `float`
            Object distance in AU.

        d_star : `float`
            Star diameter, in km.


    Warning
    -------
    Input data must be one of the 4 options below:

    1) Input data from file with time and flux
        `file (str)`: a file with the time and flux. A third column with the error in
        flux can also be given.

        `usecols (int, tuple, array)`: Which columns to read, with the first being the
        time, the seconds the flux and third the flux error (optional).

        **Example:**

        >>> LightCurve(name, file, exptime) # dflux can also be given


    2) Input data when file is not given:
        `time`: time must be a list of times, in seconds from tref, or Julian Date, or
        a Time object.

        `flux`: flux must be a list of fluxes. It must have the same lenght as time.

        `dflux`: if file not given, dflux must be a list of fluxes errors. It must
        have the same lenght as time. (not required)

        **Example:**

        >>> LightCurve(name, flux, time, exptime) # dflux can also be given


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

    """

    @deprecated_alias(lambda_0='central_bandpass', delta_lambda='delta_bandpass')  # remove this line for v1.0
    def __init__(self, name='', **kwargs):

        allowed_kwargs = ['emersion', 'emersion_err', 'immersion', 'immersion_err', 'initial_time', 'end_time',
                          'file', 'time', 'flux', 'exptime', 'central_bandpass', 'delta_bandpass', 'tref', 'dflux',
                          'skiprows', 'usecols', 'dist', 'vel', 'd_star']
        input_tests.check_kwargs(kwargs, allowed_kwargs=allowed_kwargs)
        input_done = False
        self.dflux = None
        self._name = name
        self.flux = None
        self.time_model = None
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
        if not input_done:
            try:
                self.set_flux(**kwargs)
            except:
                raise ValueError('No allowed input conditions satisfied. Please refer to the tutorial.')
        self.set_filter(central_bandpass=kwargs.get('central_bandpass', 0.70),
                        delta_bandpass=kwargs.get('delta_bandpass', 0.30))
        self.dt = 0.0

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
                self.immersion = Time(value, format='jd')
            elif hasattr(self, 'tref'):
                self.immersion = self.tref + value*u.s
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
                self.emersion = Time(value, format='jd')
            elif hasattr(self, 'tref'):
                self.emersion = self.tref + value*u.s
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
                self.initial_time = Time(value, format='jd')
            elif hasattr(self, 'tref'):
                self.initial_time = self.tref + value*u.s
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
                self.end_time = Time(value, format='jd')
            elif hasattr(self, 'tref'):
                self.end_time = self.tref + value*u.s
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
            return (self._time - self.tref).sec
        except:
            raise AttributeError("'LightCurve' object has no attribute 'time'")

    def set_flux(self, **kwargs):
        """Sets the flux for the LightCurve.

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
            self.model = np.ones(len(time))
            self.flux = self.flux[order]
            self.flux_obs = self.flux
            if self.dflux is not None:
                self.dflux = self.dflux[order]
            self.initial_time = np.min(time)
            self.end_time = np.max(time)
            time_diffs = np.diff(self._time[0:]).tolist()
            self.cycle = max(set(time_diffs), key=time_diffs.count).sec
            if self.cycle < self.exptime:
                warnings.warn('Exposure time ({:0.4f} seconds) higher than Cycle time ({:0.4f} seconds)'.
                              format(self.exptime, self.cycle))

    def set_exptime(self, exptime):
        """Sets the light curve exposure time.

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
        """Sets the occultation velocity.

        Parameters
        ----------
        vel : `int`, `float`
            Velocity in km/s.
        """
        vel = u.Quantity(vel, unit=u.km/u.s)
        self.vel = np.absolute(vel.value)

    def set_dist(self, dist):
        """Sets the object distance.

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
        """Sets the star diameter.

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
        """Sets the filter bandwidth in microns.

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

    def calc_magnitude_drop(self, mag_star, mag_obj):
        """Determines the magnitude drop of the occultation.

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
        """Returns the fresnel scale.

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

        # Create a mask where the polynomial fit will be done
        if not all(self.flux):
            raise ValueError('Normalization is only possible when a LightCurve is instantiated with time and flux.')
        self.reset_flux()
        lc_flux = (self.flux - flux_min)/(flux_max-flux_min)
        if mask is None:
            preliminar_occ = self.occ_detect(maximum_duration=((self.end_time - self.initial_time).value*u.d.to('s'))/3)
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

    def occ_model(self, immersion_time, emersion_time, opacity, mask, npt_star=12,
                  time_resolution_factor=10, flux_min=0, flux_max=1):
        """Returns the modelled light curve.

        The modelled light curve takes into account the fresnel diffraction, the
        star diameter and the instrumental response.

        Parameters
        ----------
        immersion_time : `int`, `float`
            Immersion time, in seconds.

        emersion_time : `int`, `float`
            Emersion time, in seconds.

        opacity : `int`, `float`
            Opacity. Opaque = 1.0, transparent = 0.0,

        mask : `bool` array
            Mask with True values to be computed.

        npt_star : `int`, default=12
            Number of subdivisions for computing the star size effects.

        time_resolution_factor : `int`, `float`, default: 10*fresnel scale
            Steps for fresnel scale used for modelling the light curve.

        flux_min : `int`, `float`, default=0
            Bottom flux (only object).

        flux_max : `int`, `float`, default=1
            Base flux (object plus star).
        """
        from .utils import bar_fresnel

        # Computing the fresnel scale
        lamb = self.lambda_0*u.micrometer.to('km')
        dlamb = self.delta_lambda*u.micrometer.to('km')
        dist = self.dist*u.au.to('km')
        vel = np.absolute(self.vel)
        time_obs = self.time[mask]
        fresnel_scale_1 = calc_fresnel(dist, lamb-dlamb/2.0)
        fresnel_scale_2 = calc_fresnel(dist, lamb+dlamb/2.0)
        fresnel_scale = (fresnel_scale_1 + fresnel_scale_2)/2.0
        time_resolution = (np.min([fresnel_scale/vel, self.exptime]))/time_resolution_factor

        # Creating a high resolution curve to compute fresnel diffraction, stellar diameter and instrumental integration
        time_model = np.arange(time_obs.min()-5*self.exptime, time_obs.max()+5*self.exptime, time_resolution)

        # Changing X: time (s) to distances in the sky plane (km), considering the tangential velocity (vel in km/s)
        x = time_model*vel
        x01 = immersion_time*vel
        x02 = emersion_time*vel

        # Computing fresnel diffraction for the case where the star size is negligenciable
        flux_fresnel_1 = bar_fresnel(x, x01, x02, fresnel_scale_1, opacity)
        flux_fresnel_2 = bar_fresnel(x, x01, x02, fresnel_scale_2, opacity)
        flux_fresnel = (flux_fresnel_1 + flux_fresnel_2)/2.
        flux_star = flux_fresnel.copy()
        if self.d_star > 0:
            # Computing fresnel diffraction for the case where the star size is not negligenciable
            resolucao = (self.d_star/2)/npt_star
            flux_star_1 = np.zeros(len(time_model))
            flux_star_2 = np.zeros(len(time_model))
            # Computing stellar diameter only near the immersion or emersion times
            star_diam = (np.absolute(x - x01) < 3*self.d_star) + (np.absolute(x - x02) < 3*self.d_star)
            p = np.arange(-npt_star, npt_star)*resolucao
            coeff = np.sqrt(np.absolute((self.d_star/2)**2 - p**2))
            for ii in np.where(star_diam == True)[0]:
                xx = x[ii] + p
                flux1 = bar_fresnel(xx, x01, x02, fresnel_scale_1, opacity)
                flux2 = bar_fresnel(xx, x01, x02, fresnel_scale_2, opacity)
                flux_star_1[ii] = np.sum(coeff*flux1)/coeff.sum()
                flux_star_2[ii] = np.sum(coeff*flux2)/coeff.sum()
                flux_star[ii] = (flux_star_1[ii] + flux_star_2[ii])/2.
        flux_inst = np.zeros(len(time_obs))
        for i in range(len(time_obs)):
            event_model = (time_model > time_obs[i]-self.exptime/2.) & (time_model < time_obs[i]+self.exptime/2.)
            flux_inst[i] = (flux_star[event_model]).mean()
        self.model[mask] = flux_inst*(flux_max - flux_min) + flux_min
        self.time_model = time_model
        self.model_star = flux_star*(flux_max - flux_min) + flux_min
        self.model_fresnel = flux_fresnel*(flux_max - flux_min) + flux_min
        ev_model = (time_model > immersion_time) & (time_model < emersion_time)
        flux_box = np.ones(len(time_model))
        flux_box[ev_model] = (1-opacity)**2
        flux_box = flux_box*(flux_max - flux_min) + flux_min
        self.model_geometric = flux_box
        self.baseflux = flux_max
        self.bottomflux = flux_min

    def occ_lcfit(self, **kwargs):
        """Monte Carlo chi square fit for occultations lightcurve.

        Parameters
        ----------
        tmin : `int`, `float`
            Minimum time to consider in the fit procedure, in seconds.

        tmax : `int`, `float`
            Maximum time to consider in the fit procedure, in seconds.

        flux_min : `int`, `float`, default=0
            Bottom flux (only object).

        flux_max :`int`, `float`, default=1
            Base flux (object plus star).

        immersion_time : `int`, `float`
            Initial guess for immersion time, in seconds.

        emersion_time : `int`, `float`
            Initial guess for emersion time, in seconds.

        opacity : `int`, `float`, default=1
            Initial guess for opacity. Opaque = 1, Transparent = 0.

        delta_t : `int`, `float`
            Interval to fit immersion or emersion time.

        dopacity : `int`, `float`, default=0
            Interval to fit opacity.

        sigma : `int`, `float`, `array`, 'auto'
            Fluxes errors. If None it will use the `self.dflux`. If 'auto' it
            will calculate using the region outside the event.

        loop : `int`, default=10000
            Number of tests to be done.

        sigma_result : `int`, `float`
            Sigma value to be considered as result.

        Returns
        -------
        chi2 : `sora.extra.ChiSquare`
            ChiSquare object.
        """
        from sora.config.visuals import progressbar
        from sora.extra import ChiSquare

        allowed_kwargs = ['tmin', 'tmax', 'flux_min', 'flux_max', 'immersion_time', 'emersion_time', 'opacity',
                          'delta_t', 'dopacity', 'sigma', 'loop', 'sigma_result']
        input_tests.check_kwargs(kwargs, allowed_kwargs=allowed_kwargs)

        if not hasattr(self, 'flux'):
            raise ValueError('Fit curve is only possible when a LightCurve is instantiated with time and flux.')

        preliminar_occ = self.occ_detect()

        delta_t = 2*self.cycle
        loop = kwargs.get('loop', 10000)
        tmax = self.time.max()
        tmin = self.time.min()
        immersion_time = tmin - self.exptime
        do_immersion = False
        emersion_time = tmax + self.exptime
        do_emersion = False
        opacity = kwargs.get('opacity', 1.0)
        delta_opacity = kwargs.get('dopacity', 0.0)
        do_opacity = 'dopacity' in kwargs
        if ('immersion_time' not in kwargs) and ('emersion_time' not in kwargs):
            immersion_time = preliminar_occ['immersion_time']
            do_immersion = True
            emersion_time = preliminar_occ['emersion_time']
            do_emersion = True
            delta_t = 5*preliminar_occ['time_err']
            tmax = emersion_time+2*preliminar_occ['occultation_duration']
            tmin = immersion_time-2*preliminar_occ['occultation_duration']
            if 2*preliminar_occ['occultation_duration'] < 10*self.cycle:
                tmax = emersion_time + 10*self.cycle
                tmin = immersion_time - 10*self.cycle
        tmax = kwargs.get('tmax', tmax)
        tmin = kwargs.get('tmin', tmin)
        delta_t = kwargs.get('delta_t', delta_t)
        if 'immersion_time' in kwargs:
            immersion_time = kwargs['immersion_time']
            do_immersion = True
        t_i = immersion_time + delta_t*(2*np.random.random(loop) - 1)
        if 'emersion_time' in kwargs:
            emersion_time = kwargs['emersion_time']
            do_emersion = True
        t_e = emersion_time + delta_t*(2*np.random.random(loop) - 1)
        mask = (self.time >= tmin) & (self.time <= tmax)
        if 'sigma' not in kwargs:
            if self.dflux is not None:
                sigma = self.dflux
            else:
                sigma = 'auto'
        else:
            if type(kwargs['sigma']) in [float, int]:
                sigma = np.repeat(kwargs['sigma'], len(self.flux))
            elif kwargs['sigma'] is None:
                sigma = self.dflux
            else:
                sigma = kwargs['sigma']
        if type(sigma) is str and sigma == 'auto':
            mask_sigma = (((self.time >= tmin) & (self.time < immersion_time - self.exptime)) +
                          ((self.time > emersion_time + self.exptime) & (self.time <= tmax)))
            sigma = np.repeat(self.flux[mask_sigma].std(ddof=1), len(self.flux))
        opas = opacity + delta_opacity*(2*np.random.random(loop) - 1)
        opas[opas > 1.], opas[opas < 0.] = 1.0, 0.0
        flux_min = kwargs.get('flux_min', 1 - preliminar_occ['depth'])
        flux_max = kwargs.get('flux_max', preliminar_occ['baseline'])
        sigma_result = kwargs.get('sigma_result', 1)

        tflag = np.zeros(loop)
        tflag[t_i > t_e] = t_i[t_i > t_e]
        t_i[t_i > t_e] = t_e[t_i > t_e]
        t_e[t_i > t_e] = tflag[t_i > t_e]
        chi2 = 999999*np.ones(loop)
        for i in progressbar(range(loop), 'LightCurve fit:'):
            model_test = self.__occ_model(t_i[i], t_e[i], opas[i], mask, flux_min=flux_min, flux_max=flux_max)
            chi2[i] = np.sum(((self.flux[mask] - model_test)**2)/(sigma[mask]**2))
        kkwargs = {}
        if do_immersion:
            kkwargs['immersion'] = t_i
        if do_emersion:
            kkwargs['emersion'] = t_e
        if do_opacity:
            kkwargs['opacity'] = opas
        chisquare = ChiSquare(chi2, len(self.flux[mask]), **kkwargs)
        result_sigma = chisquare.get_nsigma(sigma=sigma_result)
        if 'immersion' in result_sigma:
            self._immersion = self.tref + result_sigma['immersion'][0]*u.s
            self.immersion_err = result_sigma['immersion'][1]
            immersion_time = result_sigma['immersion'][0]
        else:
            try:
                immersion_time = (self._immersion.jd - self.tref.jd)*u.d.to('s')
            except:
                pass
        if 'emersion' in result_sigma:
            self._emersion = self.tref + result_sigma['emersion'][0]*u.s
            self.emersion_err = result_sigma['emersion'][1]
            emersion_time = result_sigma['emersion'][0]
        else:
            try:
                emersion_time = (self._emersion.jd - self.tref.jd)*u.d.to('s')
            except:
                pass
        if 'opacity' in result_sigma:
            opacity = result_sigma['opacity'][0]
        # Run occ_model() to save best parameters in the Object.
        self.occ_model(immersion_time, emersion_time, opacity, np.repeat(True, len(self.flux)),
                       flux_min=flux_min, flux_max=flux_max)
        self.lc_sigma = sigma
        self.chisquare = chisquare
        self.opacity = opacity
        return chisquare

    def plot_lc(self, ax=None):
        """ Plots the light curve
        """
        import matplotlib.pyplot as plt

        if not any(self.flux):
            raise ValueError('Plotting the light curve is only possible when the '
                             'Object LightCurve is instantiated with time and flux')
        ax = ax or plt.gca()
        ax.plot(self.time, self.flux, 'k.-', label='Obs.', zorder=0)
        if any(self.model):
            ax.plot(self.time, self.model, 'r-', label='Model', zorder=2)
            ax.scatter(self.time, self.model, s=50, facecolors='none', edgecolors='r', zorder=3)
        ax.set_xlabel('Time [seconds]', fontsize=20)
        ax.set_ylabel('Relative Flux', fontsize=20)
        ax.legend()

    def plot_model(self, ax=None):
        """ Plots the modelled light curve
        """
        import matplotlib.pyplot as plt

        if not all(self.time_model):
            raise ValueError('Plotting the model light curve is only possible after the model '
                             '[LightCurve.occ_model()] or the fit [LightCurve.occ_lcfit()]')
        ax = ax or plt.gca()
        ax.plot(self.time_model, self.model_geometric, 'c-', label='Geometric', zorder=1)
        ax.plot(self.time_model, self.model_fresnel, 'b-', label='Fresnel', zorder=1)
        ax.plot(self.time_model, self.model_star, 'g-', label='Star diam.', zorder=1)
        ax.set_xlabel('Time [seconds]', fontsize=20)
        ax.set_ylabel('Relative Flux', fontsize=20)
        ax.legend()

    def to_log(self, namefile=None):
        """Saves the light curve log to a file.

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

    def to_file(self, namefile=None):
        """Saves the light curve to a file.

        Parameters
        ----------
        namefile : `str`
            Filename to save the data.
        """
        # Observational data
        if namefile is None:
            folder = ''
            file = self.name.replace(' ', '_')+'.dat'
        else:
            folder = os.path.dirname(namefile)
            file = os.path.basename(namefile)
        data = np.array([(self.time*u.s + self.tref).jd, self.time, self.flux, self.model, self.flux-self.model])
        colunm_names = ['Time JD', 'Time relative to {} UTC in seconds'.format(self.tref.iso),
                        'Observational Flux', 'Modelled Flux', 'Residual O-C']
        np.savetxt(os.path.join(folder, file), data.T, fmt='%11.8f')
        f = open(os.path.join(folder, file) + '.label', 'w')
        for i, name in enumerate(colunm_names):
            f.write('Column {}: {}\n'.format(i+1, name))
        f.close()
        # Complete Model
        if all(self.time_model):
            data_model = np.array([(self.time_model*u.s + self.tref).jd, self.time_model, self.model_geometric,
                                   self.model_fresnel, self.model_star])
            colunm_names_model = ['Model time JD', 'Model time relative to {} UTC in seconds'.format(self.tref.iso),
                                  'Geometric Model', 'Model with Fresnel diffraction', 'Model with star diameter']
            np.savetxt(os.path.join(folder, 'model_'+file), data_model.T, fmt='%11.8f')
            f = open(os.path.join(folder, 'model_'+file)+'.label', 'w')
            for i, name in enumerate(colunm_names_model):
                f.write('Column {}: {}\n'.format(i+1, name))
            f.close()

    def occ_detect(self, maximum_duration=None, dur_step=None, snr_limit=None,
                   n_detections=None, plot=False):
        """Detects automatically the occultation event in the light curve.

        Detects a 'square well' shaped transit. All parameters are optional.

        Parameters
        ----------
        maximum_duration : `float`, default: light curve time span
            Maximum duration of the occultation event.

        dur_step : `float`, default: 1/2 of sampling rate
            Step size to sweep occultation duration event.

        snr_limit : `float`, default=None
            Minimum occultation SNR.

        n_detections : `int`, default=1
            Number of detections regardless the SNR. `n_detections` is
            superseded by `snr_limit`.

        plot : `bool`
            True if output plots are desired.


        Returns
        -------
        OrderedDict : `dict`
            An ordered dictionary of :attr:`name`::attr:`value` pairs for each
            parameter.

        Examples
        --------
        >>> lc = LightCurve(time=time, flux=flux, exptime=0.0, name='lc_example')
        >>> params = lc.occ_detect()
        >>> params
        {'rank': 1,
        'occultation_duration': 40.1384063065052,
        'central_time': 7916.773870512843,
        'immersion_time': 7896.7046673595905,
        'emersion_time': 7936.843073666096,
        'time_err': 0.05011036992073059,
        'depth': 0.8663887801707082,
        'depth_err': 0.10986223384336465,
        'baseline': 0.9110181732552853,
        'baseline_err': 0.19045768512595365,
        'snr': 7.886138392251848,
        'occ_mask': array([False, False, False, ..., False, False, False])}
        """
        from .occdetect import occ_detect
        occ = occ_detect(self.flux, self.dflux, self.time, self.cycle, maximum_duration=maximum_duration,
                         dur_step=dur_step, snr_limit=snr_limit, n_detections=n_detections, plot=plot)
        return occ

    def __occ_model(self, immersion_time, emersion_time, opacity, mask, npt_star=12,
                    time_resolution_factor=10, flux_min=0.0, flux_max=1.0):
        """Private function. Returns the modelled light curve.

        Returns the modelled light curve considering fresnel diffraction, star
        diameter and instrumental response, intended for fitting inside the
        `self.occ_lcfit()`.

        Parameters
        ----------
        immersion_time : `int`, `float`
            Immersion time, in seconds.

        emersion_time : `int`, `float`
            Emersion time, in seconds.

        opacity `int`, `float`
            Opacity. Opaque = 1, Transparent = 0.

        mask : `bool` array
            Mask with True values to be computed.

        npt_star : `int`, default=12
            Number of subdivisions for computing the star size effects.

        time_resolution_factor : `int`, `float`, default=10*fresnel scale
            Steps for fresnel scale used for modelling the light curve.

        flux_min : `int`, `float`, default=0
            Bottom flux (only object).

        flux_max : `int`, `float`, default=1
            Base flux (object plus star).

        Returns
        -------
        flux_inst : array
            Modelled Instrumental light flux.
        """
        from .utils import bar_fresnel

        # Computing the fresnel scale
        lamb = self.lambda_0*u.micrometer.to('km')
        dlamb = self.delta_lambda*u.micrometer.to('km')
        dist = self.dist*u.au.to('km')
        vel = np.absolute(self.vel)
        time_obs = self.time[mask]
        fresnel_scale_1 = calc_fresnel(dist, lamb-dlamb/2.0)
        fresnel_scale_2 = calc_fresnel(dist, lamb+dlamb/2.0)
        fresnel_scale = (fresnel_scale_1 + fresnel_scale_2)/2.0
        time_resolution = (np.min([fresnel_scale/vel, self.exptime]))/time_resolution_factor
        self.model_resolution = time_resolution

        # Creating a high resolution curve to compute fresnel diffraction, stellar diameter and instrumental integration
        time_model = np.arange(time_obs.min()-5*self.exptime, time_obs.max()+5*self.exptime, time_resolution)

        # Changing X: time (s) to distances in the sky plane (km), considering the tangential velocity (vel in km/s)
        x = time_model*vel
        x01 = immersion_time*vel
        x02 = emersion_time*vel

        # Computing fresnel diffraction for the case where the star size is negligenciable
        flux_fresnel_1 = bar_fresnel(x, x01, x02, fresnel_scale_1, opacity)
        flux_fresnel_2 = bar_fresnel(x, x01, x02, fresnel_scale_2, opacity)
        flux_fresnel = (flux_fresnel_1 + flux_fresnel_2)/2.
        flux_star = flux_fresnel.copy()
        if self.d_star > 0:
            # Computing fresnel diffraction for the case where the star size is not negligenciable
            resolucao = (self.d_star/2)/npt_star
            flux_star_1 = np.zeros(len(time_model))
            flux_star_2 = np.zeros(len(time_model))
            # Computing stellar diameter only near the immersion or emersion times
            star_diam = (np.absolute(x - x01) < 3*self.d_star) + (np.absolute(x - x02) < 3*self.d_star)
            p = np.arange(-npt_star, npt_star)*resolucao
            coeff = np.sqrt(np.absolute((self.d_star/2)**2 - p**2))
            for ii in np.where(star_diam == True)[0]:
                xx = x[ii] + p
                flux1 = bar_fresnel(xx, x01, x02, fresnel_scale_1, opacity)
                flux2 = bar_fresnel(xx, x01, x02, fresnel_scale_2, opacity)
                flux_star_1[ii] = np.sum(coeff*flux1)/coeff.sum()
                flux_star_2[ii] = np.sum(coeff*flux2)/coeff.sum()
                flux_star[ii] = (flux_star_1[ii] + flux_star_2[ii])/2.
        flux_inst = np.zeros(len(time_obs))
        for i in range(len(time_obs)):
            event_model = (time_model > time_obs[i]-self.exptime/2.) & (time_model < time_obs[i]+self.exptime/2.)
            flux_inst[i] = (flux_star[event_model]).mean()
        return flux_inst*(flux_max - flux_min) + flux_min

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
            output += ('Bandpass:             {:.3f} +/- {:.3f} microns\n'
                       'Object Distance:      {:.2f} AU\n'
                       'Used shadow velocity: {:.3f} km/s\n'
                       'Fresnel scale:        {:.3f} seconds or {:.2f} km\n'
                       'Stellar size effect:  {:.3f} seconds or {:.2f} km\n'.format(
                           self.lambda_0, self.delta_lambda, self.dist, self.vel,
                           self.fresnel_scale/self.vel, self.fresnel_scale,
                           self.d_star/self.vel, self.d_star)
                       )
        except:
            output += '\nThere is no occultation associated with this light curve.\n'
        try:
            output += ('Inst. response:       {:.3f} seconds or {:.2f} km\n'
                       'Dead time effect:     {:.3f} seconds or {:.2f} km\n'
                       'Model resolution:     {:.3f} seconds or {:.2f} km\n'
                       'Modelled baseflux:    {:.3f}\n'
                       'Modelled bottomflux:  {:.3f}\n'
                       'Light curve sigma:    {:.3f}\n\n'.format(
                           self.exptime, self.exptime*self.vel, self.cycle-self.exptime,
                           (self.cycle-self.exptime)*self.vel, self.model_resolution,
                           self.model_resolution*self.vel, self.baseflux, self.bottomflux,
                           self.lc_sigma.mean())
                       )
        except:
            output += '\nObject LightCurve model was not fitted.\n\n'
        try:
            output += ('Immersion time: {} UTC +/- {:.3f} seconds\n'
                       'Emersion time:  {} UTC +/- {:.3f} seconds\n\n'.format(
                           self.immersion.iso, self.immersion_err,
                           self.emersion.iso, self.emersion_err)
                       )
        except:
            output += 'Immersion and emersion times were not fitted or instantiated.\n\n'

        try:
            output += 'Monte Carlo chi square fit.\n\n' + self.chisquare.__str__() + '\n'
        except:
            pass
        return output
