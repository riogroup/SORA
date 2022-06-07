import numpy as np
import astropy.units as u
from .utils import bar_fresnel, calc_fresnel

def occ_model_fit(time, immersion_time, emersion_time, opacity, 
                  central_bandpass, delta_bandpass, distance, velocity, exptime, star_diameter,
                  npt_star=12, time_resolution_factor=10, flux_min=0, flux_max=1):
    """Returns the model of the light curve.

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

    central_bandpass : `int`, `float`, otpional, default=0.7
        The center band pass of the detector used in observation. Value in microns.

    delta_bandpass : `int`, `float`, optional, default=0.3
        The band pass width of the detector used in observation. Value in microns.

    distance : `int`, `float`:
        Object distance in AU.
    
    velocity : `int`, `float`
        Velocity in km/s.

    exptime : `int`, `float`
        The exposure time of the observation, in seconds.

    star_diameter : `float`
        Star diameter, in km.        

    npt_star : `int`, default=12
        Number of subdivisions for computing the star size effects.

    time_resolution_factor : `int`, `float`, default: 10*fresnel scale
        Steps for fresnel scale used for modelling the light curve.

    flux_min : `int`, `float`, default=0
        Bottom flux (only object).

    flux_max : `int`, `float`, default=1
        Base flux (object plus star).
    """

    # Computing the fresnel scale
    lamb = central_bandpass*u.micrometer.to('km')
    dlamb = delta_bandpass*u.micrometer.to('km')
    dist = distance*u.au.to('km')
    vel = np.absolute(velocity)
    time_obs = time
    fresnel_scale_1 = calc_fresnel(dist, lamb-dlamb/2.0)
    fresnel_scale_2 = calc_fresnel(dist, lamb+dlamb/2.0)
    fresnel_scale = (fresnel_scale_1 + fresnel_scale_2)/2.0
    time_resolution = (np.min([fresnel_scale/vel, exptime]))/time_resolution_factor

    # Creating a high resolution curve to compute fresnel diffraction, stellar diameter and instrumental integration
    time_model = np.arange(time_obs.min()-5*exptime, time_obs.max()+5*exptime, time_resolution)

    # Changing X: time (s) to distances in the sky plane (km), considering the tangential velocity (vel in km/s)
    x = time_model*vel
    x01 = immersion_time*vel
    x02 = emersion_time*vel

    # Computing fresnel diffraction for the case where the star size is negligenciable
    flux_fresnel_1 = bar_fresnel(x, x01, x02, fresnel_scale_1, opacity)
    flux_fresnel_2 = bar_fresnel(x, x01, x02, fresnel_scale_2, opacity)
    flux_fresnel = (flux_fresnel_1 + flux_fresnel_2)/2.
    flux_star = flux_fresnel.copy()
    if star_diameter > 0:
        # Computing fresnel diffraction for the case where the star size is not negligenciable
        resolucao = (star_diameter/2)/npt_star
        flux_star_1 = np.zeros(len(time_model))
        flux_star_2 = np.zeros(len(time_model))
        # Computing stellar diameter only near the immersion or emersion times
        star_diam = (np.absolute(x - x01) < 3*star_diameter) + (np.absolute(x - x02) < 3*star_diameter)
        p = np.arange(-npt_star, npt_star)*resolucao
        coeff = np.sqrt(np.absolute((star_diameter/2)**2 - p**2))
        for ii in np.where(star_diam == True)[0]:
            xx = x[ii] + p
            flux1 = bar_fresnel(xx, x01, x02, fresnel_scale_1, opacity)
            flux2 = bar_fresnel(xx, x01, x02, fresnel_scale_2, opacity)
            flux_star_1[ii] = np.sum(coeff*flux1)/coeff.sum()
            flux_star_2[ii] = np.sum(coeff*flux2)/coeff.sum()
            flux_star[ii] = (flux_star_1[ii] + flux_star_2[ii])/2.
    flux_inst = np.zeros(len(time_obs))
    for i in range(len(time_obs)):
        event_model = (time_model > time_obs[i]-exptime/2.) & (time_model < time_obs[i]+exptime/2.)
        flux_inst[i] = (flux_star[event_model]).mean()
    return flux_inst*(flux_max - flux_min) + flux_min



def occ_model_fitError(parameters, time, flux, dflux, flux_min, flux_max,
                       central_bandpass, delta_bandpass, distance, velocity, exptime, star_diameter, 
                       time_resolution_factor, npt_star):
    '''Returns the residuals when using occ_model_fit 
    
    Parameters
    ----------
    parameters : `object`
        `Parameters` object from `Stats` module.

    time: `float` array
        Time variable.

    flux: `float` array
        Flux variable

    dflux: `float` array
        Flux uncertainty variable

    mask : `bool` array
        Mask with True values to be computed.

    central_bandpass : `int`, `float`, otpional, default=0.7
        The center band pass of the detector used in observation. Value in microns.

    delta_bandpass : `int`, `float`, optional, default=0.3
        The band pass width of the detector used in observation. Value in microns.

    distance : `int`, `float`:
        Object distance in AU.
    
    velocity : `int`, `float`
        Velocity in km/s.

    exptime : `int`, `float`
        The exposure time of the observation, in seconds.

    star_diameter : `float`
        Star diameter, in km.        

    npt_star : `int`, default=12
        Number of subdivisions for computing the star size effects.

    time_resolution_factor : `int`, `float`, default: 10*fresnel scale
        Steps for fresnel scale used for modelling the light curve.

    flux_min : `int`, `float`, default=0
        Bottom flux (only object).

    flux_max : `int`, `float`, default=1
        Base flux (object plus star).
    """
    '''
    v = parameters.valuesdict()
    model = occ_model_fit(time, v['immersion_time'], v['emersion_time'], v['opacity'],  
                          central_bandpass, delta_bandpass, distance, velocity, exptime, star_diameter,
                          npt_star=npt_star, time_resolution_factor=time_resolution_factor, flux_min=flux_min, flux_max=flux_max)
    return (flux - model)**2 / dflux**2