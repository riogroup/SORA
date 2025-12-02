import warnings

import numpy as np

warnings.simplefilter('always', UserWarning)

class DetectionLimits:
    """ Class to compute detection limits for ring-like features in a stellar occultation light curve.

    This class provides estimates of the minimum detectable opacity, optical depth, and equivalent
    width of a feature based on the noise level in the flux. If available, per-point uncertainties
    (`lc.dflux`) are used in a weighted noise estimation; otherwise, a simple standard deviation is applied.

    The computation masks out points inside the occultation interval if immersion/emersion times are defined,
    using ±1 exposure time as margin. Otherwise, all flux data is used.

    Parameters
    ----------
    lightcurve : LightCurve
        The light curve object associated with the detection limit analysis.

    Attributes
    ----------
    flux : array-like
        Flux values outside the occultation region used to estimate the detection limits.

    flux_err : array-like or None
        Flux uncertainties corresponding to `flux`, or None if unavailable or inconsistent.

    Notes
    -----
    The detection limits computed here are projected in the plane of the sky and do not include
    corrections for diffraction or projection effects (e.g., inclination of the ring plane).
    See `sora.occultation.chord.DetectionLimits` for corrected estimates.

    Examples
    --------
    >>> from sora.lightcurve import LightCurve
    >>> lc = LightCurve(file='lightcurve.dat', usecols=[0,1,2], exptime=0.1, tref='2025-01-01')
    >>> lc.detection_limits.apparent_opacity()
    0.027
    >>> lc.detection_limits.apparent_optical_depth(sigma=3)
    0.080
    >>> print(lc.detection_limits)
    """
    
    def __init__(self, lightcurve):
        """ Initialize the detection limit computation by masking data and preparing flux arrays.

        If immersion/emersion/tref are defined in the light curve, the method masks out the
        data during the occultation (±1 exposure). Otherwise, the full flux array is used.

        Parameters
        ----------
        lightcurve : LightCurve
            The input light curve object containing flux and timing data.
        """

        self.lc = lightcurve
        self.flux_err = getattr(self.lc, 'dflux', None)
        self.flux = getattr(self.lc, 'flux', None)
        immersion = getattr(self.lc, 'immersion', None)
        emersion = getattr(self.lc, 'emersion', None)
        tref = getattr(self.lc, 'tref', None)

        if self.lc.flux is None or self.lc.time is None:
            return 
        
        # Define mask to exclude data inside the occultation (with margin of 1 exposure time)
        if immersion and emersion and tref:
            imm = (immersion - tref).sec - self.lc.exptime
            eme = (emersion - tref).sec + self.lc.exptime
            mask = (self.lc.time < imm) + (self.lc.time > eme)
            self.flux = np.array(self.lc.flux)[mask]

            if self.flux_err is not None:
                self.flux_err = np.array(self.flux_err)[mask]
        else:
            self.flux = np.array(self.lc.flux)
            if self.flux_err is not None:
                self.flux_err = np.array(self.flux_err)

        if self.flux_err is not None and len(self.flux_err) != len(self.flux):
            self.flux_err = None

    def _estimate_noise(self):
        """ Estimate the 1-sigma noise level in the out-of-occultation flux.

        If flux uncertainties are available and consistent in size with the flux,
        a weighted standard deviation is used. Otherwise, the sample standard deviation
        of the flux is returned.

        Returns
        -------
        float or None
            Estimated 1-sigma noise in the normalized flux, or None if data is invalid.
        """
        
        from astropy.stats import sigma_clip

        if self.flux is None or len(self.flux) == 0:
            return None

        # Remove deep occultations automatically (5σ clipping) - first estimative
        clipped = sigma_clip(self.flux, sigma=5, maxiters=5)
        clean_flux = clipped.data[~clipped.mask]

        if len(clean_flux) < 10:
            clean_flux = self.flux

        if self.flux_err is not None:
            clean_err = self.flux_err[~clipped.mask] if len(self.flux_err) == len(self.flux) else self.flux_err
            weights = 1 / clean_err**2
            wmean = np.average(clean_flux, weights=weights)
            wvar = np.average((clean_flux - wmean)**2, weights=weights)
            return np.sqrt(wvar)

        return np.std(clean_flux, ddof=1)

    def apparent_opacity(self, sigma=1):
        """  Compute the apparent opacity detection limit.

        The apparent opacity corresponds to the minimum flux drop distinguishable from noise
        at a given sigma level.

        Parameters
        ----------
        sigma : float, optional
            The sigma level for the detection threshold (default is 1-sigma).

        Returns
        -------
        float or None
            Apparent opacity value or None if flux data is missing.
        """
        noise = self._estimate_noise()
        if noise is None:
            return None
        return noise * sigma

    
    def apparent_optical_depth(self, sigma=1):
        """ Compute the apparent optical depth detection limit.

        Derived from the apparent opacity using the relation:
        τ = -ln(1 - opacity)

        Parameters
        ----------
        sigma : float, optional
            The sigma level for the detection threshold (default is 1-sigma).

        Returns
        -------
        float or None
            Apparent optical depth value or None if flux data is missing.
        """
        noise = self._estimate_noise()
        if noise is None:
            return None
        return -np.log(1 - noise*sigma)
    
    def apparent_equivalent_width(self, sigma=1):
        """ Estimate the apparent equivalent width detection limit (in km).

        Computed as the product of the flux noise, the exposure time, and the projected velocity:
        EW = sigma * flux_noise * velocity * exptime

        Parameters
        ----------
        sigma : float, optional
            The sigma level for the detection threshold (default is 1-sigma).

        Returns
        -------
        float or None
            Apparent equivalent width in km or None if flux data is missing.
        """
        noise = self._estimate_noise()
        if noise is None:
            return None
        return noise * sigma * self.lc.vel * self.lc.exptime 
    
    def __str__(self):
        """ String summary of detection limits at 3-sigma.

        Returns
        -------
        str
            Formatted output with apparent opacity, optical depth, and equivalent width.
        """
        if self.lc.flux is not None and len(self.lc.flux) > 0:
            output = '\nDetection limits (3-sigma):\n' 
            output += ('    Apparent opacity:             {:.3f}\n'
                       '    Apparent optical depth:       {:.3f}\n'
                       '    Apparent equivalent width:    {:.3f} km\n\n'.format(
                           self.apparent_opacity(sigma=3),
                           self.apparent_optical_depth(sigma=3),
                           self.apparent_equivalent_width(sigma=3)
                           )
                        )
            return output
        return '\nNo flux data available to compute detection limits.\n'