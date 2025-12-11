# sora/lightcurve/model.py
"""
Light-curve models for stellar occultations.

Features
--------
- SquareWellModel : Fresnel-aware square-well occultation model with optional spectral integration.
- CompositeModel : multiplicative (non-coherent) combination of independent structures.

Design goals
------------
- Preserve the physics of the original SORA `occ_model`:
  Fresnel diffraction, finite star, and instrumental integration.
- Add spectral integration (n_lambda or bandpass response curve).
- Keep API simple and backward-compatible (e.g., lc.square_well_model()).
- Modular architecture to allow future models (Lorentzian, Gaussian, etc).

Examples
--------
>>> m1 = lc.square_well_model(immersion=23474, emersion=23476, opacity=0.04)
>>> m2 = lc.square_well_model(immersion=23664.3, emersion=23710, opacity=1.0)
>>> m3 = lc.square_well_model(immersion=23911, emersion=23912.3, opacity=0.7)

Combine multiple components (non-coherent product of intensities):
>>> multi = CompositeModel()
>>> multi.add_component("main_body", m1)
>>> multi.add_component("ring1", m2)
>>> multi.add_component("ring2", m3)
>>> flux = multi()
>>> lc.plot(model=flux)
"""

from __future__ import annotations
import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt
from typing import Dict, Tuple, Optional
from scipy.special import fresnel
from astropy.time import Time
import os



from .utils import calc_fresnel, bar_fresnel

__all__ = [
    "BaseModel",
    "SquareWellModel",
    "CompositeModel",
    "DoubleSquareWellModel",
    "attach_to_lightcurve_class",
    "occ_model", 
    "occ_model_double"
]

# Base class
class BaseModel:
    """Base interface for light-curve models."""

    def __init__(self, name: str = "BaseModel", params: Optional[Dict] = None):
        self.name = name
        self._params: Dict = params or {}
        self.lightcurve = None

    def __call__(self) -> np.ndarray:
        """Shortcut for self.compute()."""
        return self.compute()
    
    @property
    def params(self):
        return self._params
    
    @params.setter
    def params(self, value: Dict):
        if not isinstance(value, dict):
            raise TypeError("params must be a dict.")
        self._params = value

# Fresnel square-well model with spectral integration
class SquareWellModel(BaseModel):
    """
    Square-well occultation model.

    Includes Fresnel diffraction, finite stellar diameter, instrumental integration,
    and optional spectral integration across the bandpass.

    Parameters
    ----------
    immersion, emersion : float
        Immersion and emersion times (s, relative to tref or arbitrary zero).
    opacity : float
        Opacity (1 = opaque, 0 = transparent).
    distance : float, optional
        Object-observer distance in AU.
    vel : float, optional
        Tangential vel in km/s.
    exptime : float, optional
        Exposure time in seconds.
    d_star : float, optional
        Stellar diameter (km). If 0, star is treated as a point source.
    lambda_0, delta_lambda : float, optional
        Central wavelength and bandpass width (μm).
    npt_star : int, optional
        Number of subdivisions for integrating the stellar disk.
    time_resolution_factor : float, optional
        Time sampling factor relative to Fresnel scale (default 10).
    flux_min, flux_max : float, optional
        Flux scaling limits (default 0-1).
    n_lambda : int, optional
        Number of wavelength samples across the bandpass (default 2).
        n_lambda = 1 → monochromatic (λ₀),
        n_lambda = 2 → average of λ₀ ± Δλ/2 (SORA classic),
        n_lambda > 2 → uniform integration across band.
    response : (μm, Tλ) array pair, optional
        Transmission curve for realistic spectral integration.
    lightcurve : LightCurve, optional
        If provided, inherits distance, vel, exptime, d_star, λ₀, Δλ from it.

    Notes
    -----
    This reproduces the physics of the classic `occ_model()` but in modular form.
    """

    def __init__(
            self,
            immersion: float,
            emersion: float,
            opacity: float,
            distance: Optional[float] = None,
            vel: Optional[float] = None,
            exptime: Optional[float] = None,
            d_star: Optional[float] = None,
            lambda_0: Optional[float] = None,
            delta_lambda: Optional[float] = None,
            lightcurve=None,
            npt_star: int = 12,
            time_resolution_factor: float = 10.0,
            flux_min: float = 0.0,
            flux_max: float = 1.0,
            n_lambda: int = 2,
            response: Optional[Tuple[np.ndarray, np.ndarray]] = None,
        ):
        if lightcurve is not None:
            distance     = getattr(lightcurve, "dist", distance)
            vel     = getattr(lightcurve, "vel", vel)
            exptime      = getattr(lightcurve, "exptime", exptime)
            d_star       = getattr(lightcurve, "d_star", d_star)
            lambda_0     = getattr(lightcurve, "lambda_0", lambda_0)
            delta_lambda = getattr(lightcurve, "delta_lambda", delta_lambda)
            response     = getattr(lightcurve, "response", response)
            n_lambda     = getattr(lightcurve, "n_lambda", n_lambda)


        if distance is None or vel is None or exptime is None:
            raise ValueError(
                "Missing required parameters: distance, vel, or exptime.\n"
                "Provide them explicitly or pass a LightCurve via 'lightcurve'."
            )

        super().__init__(name="SquareWell")
        self.lightcurve = lightcurve
        self.params = dict(
            immersion=float(immersion),
            emersion=float(emersion),
            opacity=float(opacity),
            distance=float(distance),
            vel=float(vel),
            exptime=float(exptime),
            d_star=float(d_star),
            lambda_0=float(lambda_0 or 0.58),
            delta_lambda=float(delta_lambda or 0.40),
            npt_star=int(npt_star),
            time_resolution_factor=float(time_resolution_factor),
            flux_min=float(flux_min),
            flux_max=float(flux_max),
            n_lambda=int(n_lambda),
            response=response,
        )

        self.model_flux = None
        self.time_model = None
        self.model_fresnel = None
        self.model_star = None
        self.model_geometric = None

        self.immersion_err = None
        self.emersion_err  = None
        self.opacity_err   = None

    @property
    def immersion(self):
        """
        Immersion as absolute Time if lightcurve.tref exists,
        otherwise returns the relative immersion (float, s).
        """
        t_rel = self.params.get("immersion")
        if t_rel is None:
            return None
        tref = getattr(self.lightcurve, "tref", None)
        return tref + t_rel * u.s

    @property
    def emersion(self):
        """
        Emersion as absolute Time if lightcurve.tref exists,
        otherwise returns the relative emersion (float, s).
        """
        t_rel = self.params.get("emersion")
        if t_rel is None:
            return None
        tref = getattr(self.lightcurve, "tref", None)
        return tref + t_rel * u.s
    
    @property
    def opacity(self):
        return self.params.get("opacity")
    
    # Main compute method
    def compute(self, time: Optional[np.ndarray] = None) -> np.ndarray:

        p = self.params

        if time is None:
            time = np.asarray(self.lightcurve.time, dtype=float)
        else:
            time = np.asarray(time, dtype=float)
 
        lam_array, weights = _lambda_weights(p["lambda_0"], 
                                             p["delta_lambda"], 
                                             p["n_lambda"],
                                             p["response"])

        dist_km = p["distance"] * u.au.to("km")
        vel = abs(p["vel"])
        lamb_center = p["lambda_0"] * u.micrometer.to("km")
        fresnel_center = calc_fresnel(dist_km, lamb_center)
        time_resolution = min(fresnel_center / vel, p["exptime"]) / p["time_resolution_factor"]
        time_model = np.arange(time.min() - 5*p["exptime"], time.max() + 5*p["exptime"], time_resolution)
        x = time_model * vel
        x01 = p["immersion"] * vel
        x02 = p["emersion"] * vel

        flux_fresnel = np.zeros_like(time_model)

        for lam, w in zip(lam_array, weights):
            lam_km = lam * u.micrometer.to("km")
            fresnel_scale = calc_fresnel(dist_km, lam_km)
            flux_temp = bar_fresnel(x, x01, x02, fresnel_scale, p["opacity"])
            flux_fresnel += w * flux_temp
            
        flux_star = flux_fresnel.copy()

        if p["d_star"] > 0:
            res = (p["d_star"]/2) / p["npt_star"]
            pgrid = np.arange(-p["npt_star"], p["npt_star"]) * res
            coeff = np.sqrt(np.maximum((p["d_star"]/2)**2 - pgrid**2, 0.0))
            coeff_sum = coeff.sum() if coeff.sum() != 0 else 1.0
            coeff = coeff / coeff_sum

            mask = (np.abs(x - x01) < 3*p["d_star"]) | (np.abs(x - x02) < 3*p["d_star"])

            for ii in np.where(mask)[0]:
                xx = x[ii] + pgrid
                flux_local = np.interp(xx, x, flux_star, left=1.0, right=1.0)
                flux_star[ii] = np.sum(coeff * flux_local)

        flux_inst = np.zeros_like(time)
        for i, t in enumerate(time):
            mask_event = (time_model > (t - p["exptime"]/2.)) & (time_model < (t + p["exptime"]/2.))
            flux_inst[i] = flux_star[mask_event].mean() if mask_event.any() else flux_star[np.argmin(np.abs(time_model - t))]

        flux_inst = flux_inst * (p["flux_max"] - p["flux_min"]) + p["flux_min"]

        self.time_model = time_model
        self.model_fresnel = flux_fresnel * (p["flux_max"] - p["flux_min"]) + p["flux_min"]
        self.model_star = flux_star * (p["flux_max"] - p["flux_min"]) + p["flux_min"]
        geom = np.ones_like(time_model)
        geom[(time_model > p["immersion"]) & (time_model < p["emersion"])] = (1 - p["opacity"])**2
        self.model_geometric = geom * (p["flux_max"] - p["flux_min"]) + p["flux_min"]
        self.model_flux = flux_inst
        return flux_inst

    def plot(self, show_components=False, ax=None):
        if self.model_fresnel is None:
            self.compute()

        return _plot_model(
            ax=ax,
            lightcurve=self.lightcurve,
            flux_model=self.model_flux,
            model_geometric=self.model_geometric,
            model_fresnel=self.model_fresnel,
            model_star=self.model_star,
            time_model=self.time_model,
            show_components=show_components,
            title=f"{self.lightcurve.name}: {self.name}"
        )
    
    def to_file(self, namefile=None, overwrite=False):
        """
        Saves the modeled light curve to an ASCII file with metadata.
        """
        if namefile is None:
            namefile = f'{self.lightcurve.tref.iso[:10].replace("-", "")}_{self.lightcurve.name.replace(" ", "_").replace("-", "_")}_model.dat'

        if os.path.exists(namefile) and not overwrite:
            raise FileExistsError(f"File '{namefile}' already exists. Use overwrite=True to replace it.")
                
        header_lines = [f"SORA Model export",
                        f"Light curve name: {self.lightcurve.name}"]

        if hasattr(self, "tref"):
            header_lines.append(f"Reference time (UTC): {self.tref.isot}")
        if hasattr(self, 'immersion'):
            imm_err = getattr(self, 'immersion_err') or 0.0
            header_lines.append(f"Immersion time:     {self.immersion.iso} +/- {imm_err:.3f} s")
        if hasattr(self, 'emersion'):
            eme_err = getattr(self, 'emersion_err') or 0.0
            header_lines.append(f"Emersion time:     {self.emersion.iso} +/- {eme_err:.3f} s ")
        if hasattr(self, 'opacity'):
            opa_err = getattr(self, 'opacity_err') or 0.0
            header_lines.append(f"Opacity:     {self.opacity:.3f} +/- {opa_err:.3f}")

        header_lines.append("")
        header_lines.append("Columns: jd, sec from tref, Fresnel model, star model, geometric model")
        header = "\n".join(header_lines)

        # --- Data ---
        time_sec = self.time_model
        time_iso = Time(self.lightcurve.tref) + time_sec*u.s
        time_jd = time_iso.jd
        data = np.column_stack((time_jd, time_sec, self.model_fresnel, self.model_star, self.model_geometric))

        np.savetxt(namefile, data, header=header, fmt="%.6f")

        return namefile
    
    def __str__(self):
        string = ['-' * 79]
        string.append(f'{self.name}')
        imm_str = f"+/- {self.immersion_err:.3f}" if self.immersion_err is not None else ''
        eme_str = f"+/- {self.emersion_err:.3f}" if self.emersion_err is not None else ''
        opa_str = f"+/- {self.opacity_err:.3f}" if self.opacity_err is not None else ''
        string.append(f"  immersion = {self.immersion} {imm_str}")
        string.append(f"  emersion  = {self.emersion} {eme_str}")
        string.append(f"  opacity   = {self.opacity} {opa_str}")
        string.append('')
        return '\n'.join(string)

# Double Square-Well model (two consecutive regions with coherent diffraction)

class DoubleSquareWellModel(BaseModel):
    """
    Fresnel-aware occultation model with TWO consecutive square-well regions.

    Includes Fresnel diffraction (coherent sum of the four edges),
    finite stellar diameter convolution, instrumental integration,
    and optional spectral integration across the bandpass.

    Parameters
    ----------
    immersion_1, emersion_1 : float
        Immersion and emersion times of the first region (s, relative to tref).
    opacity_1 : float
        Opacity of the first region (1 = opaque, 0 = transparent).
    immersion_2, emersion_2 : float
        Immersion and emersion times of the second region (s).
    opacity_2 : float
        Opacity of the second region (1 = opaque, 0 = transparent).

    distance : float, optional
        Object-observer distance in AU.
    vel : float, optional
        Tangential vel in km/s.
    exptime : float, optional
        Exposure time in seconds.
    d_star : float, optional
        Stellar diameter (km). If 0, star is treated as a point source.
    lambda_0, delta_lambda : float, optional
        Central wavelength and bandpass width (μm).
    npt_star : int, optional
        Number of subdivisions for integrating the stellar disk.
    time_resolution_factor : float, optional
        Time sampling factor relative to Fresnel scale (default 10).
    flux_min, flux_max : float, optional
        Flux scaling limits (default 0–1).
    n_lambda : int, optional
        Number of wavelength samples across the bandpass (default 2).
        n_lambda = 1 → monochromatic (λ₀),
        n_lambda = 2 → average of λ₀ ± Δλ/2 (classic),
        n_lambda > 2 → uniform integration across band.
    response : (μm, Tλ) array pair, optional
        Transmission curve for realistic spectral integration.
    lightcurve : LightCurve, optional
        If provided, inherits distance, vel, exptime, d_star, λ₀, Δλ, response from it.

    Notes
    -----
    - Physics: coherent addition of edge amplitudes for two square wells,
      then |A|^2 to get intensity. Matches the legacy SORA two-box behavior.
    - Internal high-res grid is stored in self.time_model and component curves
      in self.model_fresnel, self.model_star, self.model_geometric.
    """

    def __init__(
        self,
        immersion1: float,
        emersion1: float,
        opacity1: float,
        immersion2: float,
        emersion2: float,
        opacity2: float,
        distance: Optional[float] = None,
        vel: Optional[float] = None,
        exptime: Optional[float] = None,
        d_star: Optional[float] = None,
        lambda_0: Optional[float] = None,
        delta_lambda: Optional[float] = None,
        lightcurve=None,
        npt_star: int = 12,
        time_resolution_factor: float = 10.0,
        flux_min: float = 0.0,
        flux_max: float = 1.0,
        n_lambda: int = 2,
        response: Optional[Tuple[np.ndarray, np.ndarray]] = None,
    ):
        if lightcurve is not None:
            distance     = getattr(lightcurve, "dist", distance)
            vel     = getattr(lightcurve, "vel", vel)
            exptime      = getattr(lightcurve, "exptime", exptime)
            d_star       = getattr(lightcurve, "d_star", d_star)
            lambda_0     = getattr(lightcurve, "lambda_0", lambda_0)
            delta_lambda = getattr(lightcurve, "delta_lambda", delta_lambda)
            response     = getattr(lightcurve, "response", response)

        if distance is None or vel is None or exptime is None:
            raise ValueError(
                "Missing required parameters: distance, vel, or exptime.\n"
                "Provide them explicitly or pass a LightCurve via 'lightcurve'."
            )

        super().__init__(name="DoubleSquareWell")
        self.lightcurve = lightcurve
        self.params = dict(
            immersion1=float(immersion1),
            emersion1=float(emersion1),
            opacity1=float(opacity1),
            immersion2=float(immersion2),
            emersion2=float(emersion2),
            opacity2=float(opacity2),
            distance=float(distance),
            vel=float(vel),
            exptime=float(exptime),
            d_star=float(d_star or 0.0),
            lambda_0=float(lambda_0 or 0.58),
            delta_lambda=float(delta_lambda or 0.40),
            npt_star=int(npt_star),
            time_resolution_factor=float(time_resolution_factor),
            flux_min=float(flux_min),
            flux_max=float(flux_max),
            n_lambda=int(n_lambda),
            response=response,
        )

        self.model_flux = None
        self.time_model = None
        self.model_fresnel = None
        self.model_star = None
        self.model_geometric = None

        self.immersion1_err = None
        self.emersion1_err  = None
        self.opacity1_err   = None
        self.immersion2_err = None
        self.emersion2_err  = None
        self.opacity2_err   = None


    @property
    def immersion1(self):
        """
        Immersion as absolute Time if lightcurve.tref exists,
        otherwise returns the relative immersion (float, s).
        """
        t_rel = self.params.get("immersion1")
        if t_rel is None:
            return None
        tref = getattr(self.lightcurve, "tref", None)
        return tref + t_rel * u.s

    @property
    def emersion1(self):
        """
        Emersion as absolute Time if lightcurve.tref exists,
        otherwise returns the relative emersion (float, s).
        """
        t_rel = self.params.get("emersion1")
        if t_rel is None:
            return None
        tref = getattr(self.lightcurve, "tref", None)
        return tref + t_rel * u.s
    
    @property
    def opacity1(self):
        return self.params.get("opacity1")
    
    @property
    def immersion2(self):
        """
        Immersion as absolute Time if lightcurve.tref exists,
        otherwise returns the relative immersion (float, s).
        """
        t_rel = self.params.get("immersion2")
        if t_rel is None:
            return None
        tref = getattr(self.lightcurve, "tref", None)
        return tref + t_rel * u.s

    @property
    def emersion2(self):
        """
        Emersion as absolute Time if lightcurve.tref exists,
        otherwise returns the relative emersion (float, s).
        """
        t_rel = self.params.get("emersion2")
        if t_rel is None:
            return None
        tref = getattr(self.lightcurve, "tref", None)
        return tref + t_rel * u.s
    
    @property
    def opacity2(self):
        return self.params.get("opacity2")


    def compute(self, time: Optional[np.ndarray] = None) -> np.ndarray:
        """
        Compute model flux using either the LightCurve time array or a custom one.

        Parameters
        ----------
        time : array-like, optional
            Custom time array (in seconds relative to tref).
            If not given, uses self.lightcurve.time by default.

        Returns
        -------
        flux_inst : ndarray
            Modeled relative flux values evaluated at the requested time grid.
        """

        if self.lightcurve is None and time is None:
            raise ValueError("This model needs a LightCurve or a provided 'time' array.")

        p = self.params

        if time is None:
            if not hasattr(self.lightcurve, "time"):
                raise ValueError("Linked LightCurve has no 'time' array.")
            time = np.asarray(self.lightcurve.time, dtype=float)
        else:
            time = np.asarray(time, dtype=float)

        lam_array, weights = _lambda_weights(p["lambda_0"], 
                                             p["delta_lambda"], 
                                             p["n_lambda"], 
                                             p["response"])

        dist_km = p["distance"] * u.au.to("km")
        vel = abs(p["vel"])
        lamb_center = p["lambda_0"] * u.micrometer.to("km")
        fresnel_center = calc_fresnel(dist_km, lamb_center)
        time_resolution = min(fresnel_center / vel, p["exptime"]) / p["time_resolution_factor"]
        tmin = time.min() - 5.0 * p["exptime"]
        tmax = time.max() + 5.0 * p["exptime"]
        time_model = np.arange(tmin, tmax, time_resolution)

        x = time_model * vel
        x01 = p["immersion1"] * vel
        x02 = p["emersion1"] * vel
        x03 = p["immersion2"] * vel
        x04 = p["emersion2"] * vel

        flux_fresnel = np.zeros_like(time_model)

        for lam, w in zip(lam_array, weights):
            lam_km = lam * u.micrometer.to("km")
            fr = calc_fresnel(dist_km, lam_km)

            u1 = (x - x01) / fr
            u2 = (x - x02) / fr
            u3 = (x - x03) / fr
            u4 = (x - x04) / fr

            s1, c1 = fresnel(u1)
            s2, c2 = fresnel(u2)
            s3, c3 = fresnel(u3)
            s4, c4 = fresnel(u4)

            cc1, ss1 = (c1 - c2), (s1 - s2)
            cc2, ss2 = (c3 - c4), (s3 - s4)

            r_amp = - (cc1 + ss1) * (p["opacity1"]/2.0) - (cc2 + ss2) * (p["opacity2"]/2.0)
            i_amp = + (cc1 - ss1) * (p["opacity1"]/2.0) + (cc2 - ss2) * (p["opacity2"]/2.0)

            I = (1.0 + r_amp)**2 + i_amp**2
            flux_fresnel += w * I

        flux_star = flux_fresnel.copy()

        if p["d_star"] > 0:
            npt = p["npt_star"]
            res = (p["d_star"]/2.0) / npt
            pgrid = np.arange(-npt, npt) * res
            coeff = np.sqrt(np.abs((p["d_star"]/2.0)**2 - pgrid**2))
            csum = coeff.sum() if coeff.sum() != 0 else 1.0

            edge_mask = (np.abs(x-x01) < 3*p["d_star"]) | (np.abs(x-x02) < 3*p["d_star"]) | \
                        (np.abs(x-x03) < 3*p["d_star"]) | (np.abs(x-x04) < 3*p["d_star"])
            idx = np.where(edge_mask)[0]

            for ii in idx:
                xx = x[ii] + pgrid
                acc = 0.0
                for lam, w in zip(lam_array, weights):
                    lam_km = lam * u.micrometer.to("km")
                    fr = calc_fresnel(dist_km, lam_km)

                    u1 = (xx - x01)/fr; s1, c1 = fresnel(u1)
                    u2 = (xx - x02)/fr; s2, c2 = fresnel(u2)
                    u3 = (xx - x03)/fr; s3, c3 = fresnel(u3)
                    u4 = (xx - x04)/fr; s4, c4 = fresnel(u4)

                    cc1, ss1 = (c1 - c2), (s1 - s2)
                    cc2, ss2 = (c3 - c4), (s3 - s4)

                    r_amp = - (cc1 + ss1) * (p["opacity1"]/2.0) - (cc2 + ss2) * (p["opacity2"]/2.0)
                    i_amp = + (cc1 - ss1) * (p["opacity1"]/2.0) + (cc2 - ss2) * (p["opacity2"]/2.0)

                    I = (1.0 + r_amp)**2 + i_amp**2
                    acc += w * np.sum(coeff * I) / csum

                flux_star[ii] = acc

        flux_inst = np.zeros_like(time)
        half = p["exptime"]/2.0
        for i, t in enumerate(time):
            m = (time_model > (t - half)) & (time_model < (t + half))
            flux_inst[i] = flux_star[m].mean() if m.any() else flux_star[np.argmin(np.abs(time_model - t))]

        flux_inst = flux_inst * (p["flux_max"] - p["flux_min"]) + p["flux_min"]

        self.time_model = time_model
        geom = np.ones_like(time_model)
        in1 = (time_model > p["immersion1"]) & (time_model < p["emersion1"])
        in2 = (time_model > p["immersion2"]) & (time_model < p["emersion2"])
        if np.any(in1): geom[in1] *= (1 - p["opacity1"])**2
        if np.any(in2): geom[in2] *= (1 - p["opacity2"])**2
        self.model_geometric = geom * (p["flux_max"] - p["flux_min"]) + p["flux_min"]
        self.model_fresnel   = flux_fresnel * (p["flux_max"] - p["flux_min"]) + p["flux_min"]
        self.model_star      = flux_star * (p["flux_max"] - p["flux_min"]) + p["flux_min"]
        self.model_flux      = flux_inst

        return flux_inst

    def plot(self, show_components: bool = False, ax=None):
        """Plot observed LC + geometric (always) + final model; optional Fresnel/stellar."""
        if self.model_fresnel is None:
            self.compute()
        return _plot_model(
            ax=ax,
            lightcurve=self.lightcurve,
            flux_model=self.model_flux if self.model_flux is not None else self(),
            model_geometric=self.model_geometric,
            model_fresnel=self.model_fresnel,
            model_star=self.model_star,
            time_model=self.time_model,
            show_components=show_components,
            title=f"{self.lightcurve.name}: {self.name}",
        )
    
    def to_file(self, namefile=None, overwrite=False):
        """
        Saves the modeled light curve to an ASCII file with metadata.
        """
        if namefile is None:
            namefile = f'{self.lightcurve.tref.iso[:10].replace("-", "")}_{self.lightcurve.name.replace(" ", "_").replace("-", "_")}_model.dat'

        if os.path.exists(namefile) and not overwrite:
            raise FileExistsError(f"File '{namefile}' already exists. Use overwrite=True to replace it.")
                
        header_lines = [f"SORA Model export",
                        f"Light curve name: {self.lightcurve.name}"]

        if hasattr(self, "tref"):
            header_lines.append(f"Reference time (UTC): {self.tref.isot}")
        if hasattr(self, 'immersion1'):
            imm_err = getattr(self, 'immersion1_err') or 0.0
            header_lines.append(f"Immersion1:     {self.immersion1.iso} +/- {imm_err:.3f} s")
        if hasattr(self, 'emersion1'):
            eme_err = getattr(self, 'emersion1_err') or 0.0
            header_lines.append(f"Emersion1:     {self.emersion1.iso} +/- {eme_err:.3f} s ")
        if hasattr(self, 'opacity1'):
            opa_err = getattr(self, 'opacity1_err') or 0.0
            header_lines.append(f"Opacity1:     {self.opacity1:.3f} +/- {opa_err:.3f}")
        if hasattr(self, 'immersion2'):
            imm_err = getattr(self, 'immersion2_err') or 0.0
            header_lines.append(f"Immersion2:     {self.immersion2.iso} +/- {imm_err:.3f} s")
        if hasattr(self, 'emersion1'):
            eme_err = getattr(self, 'emersion2_err') or 0.0
            header_lines.append(f"Emersion2:     {self.emersion2.iso} +/- {eme_err:.3f} s ")
        if hasattr(self, 'opacity1'):
            opa_err = getattr(self, 'opacity2_err') or 0.0
            header_lines.append(f"Opacity1:     {self.opacity2:.3f} +/- {opa_err:.3f}")

        header_lines.append("")
        header_lines.append("Columns: jd, sec from tref, Fresnel model, star model, geometric model")
        header = "\n".join(header_lines)

        # --- Data ---
        time_sec = self.time_model
        time_iso = Time(self.lightcurve.tref) + time_sec*u.s
        time_jd = time_iso.jd
        data = np.column_stack((time_jd, time_sec, self.model_fresnel, self.model_star, self.model_geometric))

        np.savetxt(namefile, data, header=header, fmt="%.6f")

        return namefile
    
    def __str__(self):
        """String summary of DoubleSquareWellModel parameters."""
        p = self.params

        lines = [f"DoubleSquareWellModel:"]

        lines.append(f"  immersion1 = {self.lightcurve.tref+ p['immersion1']*u.s}")
        lines.append(f"  emersion1  = {self.lightcurve.tref + p['emersion1']*u.s}")
        lines.append(f"  opacity1   = {p['opacity1']:.3f}")
        lines.append(f"  immersion2 = {self.lightcurve.tref+ p['immersion2']*u.s}")
        lines.append(f"  emersion2  = {self.lightcurve.tref + p['emersion2']*u.s}")
        lines.append(f"  opacity2   = {p['opacity2']:.3f}")
        return "\n".join(lines)


# Composite (non-coherent) model
class CompositeModel(BaseModel):
    """
    Multiply intensities of independent components (non-coherent).

    Generates a single combined model (geometric, Fresnel, stellar, instrumental)
    from all subcomponents.

    Example
    -------
    >>> m1 = lc.square_well_model(...)
    >>> m2 = lc.square_well_model(...)
    >>> comp = lc.composite_model()
    >>> comp.add_component("main", m1)
    >>> comp.add_component("ring1", m2)
    >>> comp.plot(show_components=True)
    """

    def __init__(self, components=None):
        super().__init__(name="CompositeModel")
        self.components = components or []
        self.lightcurve = None

        self.time_model = None
        self.model_geometric = None
        self.model_fresnel = None
        self.model_star = None
        self.model_flux = None

    def add_component(self, flag: str, model: BaseModel):
        """Add a component model to the composite."""
        if not isinstance(model, BaseModel):
            raise TypeError("Each component must be a BaseModel instance.")
        self.components.append((flag, model))
        if self.lightcurve is None:
            self.lightcurve = model.lightcurve

    def compute(self, time: Optional[np.ndarray] = None) -> np.ndarray:
        """Compute the combined model from all components.

        Parameters
        ----------
        time : array-like, optional
            Custom time array (in seconds relative to tref).
            If not given, uses the lightcurve time from the first component.

        Returns
        -------
        flux_inst : ndarray
            Combined modeled flux evaluated at the specified time array.
        """
        if not self.components:
            raise ValueError("No components added to CompositeModel.")

        lc = self.components[0][1].lightcurve
        self.lightcurve = lc

        if time is None:
            time = lc.time
        else:
            time = np.asarray(time, dtype=float)

        finest_dt = min([
            (m.params["exptime"] / m.params["time_resolution_factor"])
            for _, m in self.components
        ])
        exptime_ref = min([m.params["exptime"] for _, m in self.components])

        tmin = time.min() - 5 * exptime_ref
        tmax = time.max() + 5 * exptime_ref
        time_model = np.arange(tmin, tmax, finest_dt)

        geom_total = np.ones_like(time_model)
        fres_total = np.ones_like(time_model)
        star_total = np.ones_like(time_model)

        for _, model in self.components:
            if model.model_fresnel is None:
                model.compute()

            fres_total *= np.interp(time_model, model.time_model, model.model_fresnel)
            star_total *= np.interp(time_model, model.time_model, model.model_star)
            geom_total *= np.interp(time_model, model.time_model, model.model_geometric)

        exptime = lc.exptime
        flux_inst = np.zeros_like(time)
        for i, t in enumerate(time):
            m = (time_model > (t - exptime/2)) & (time_model < (t + exptime/2))
            flux_inst[i] = star_total[m].mean() if np.any(m) else star_total[np.argmin(np.abs(time_model - t))]

        self.time_model = time_model
        self.model_geometric = geom_total
        self.model_fresnel = fres_total
        self.model_star = star_total
        self.model_flux = flux_inst

        return flux_inst


    def plot(self, show_components=False, ax=None):
        """Plot the combined composite model and observed light curve."""
        if self.model_fresnel is None:
            self.compute()
        return _plot_model(
            ax=ax,
            lightcurve=self.lightcurve,
            flux_model=self.model_flux,
            model_geometric=self.model_geometric,
            model_fresnel=self.model_fresnel,
            model_star=self.model_star,
            time_model=self.time_model,
            show_components=show_components,
            title=f"{self.lightcurve.name}: Composite Model"
        )
    
    def __call__(self) -> np.ndarray:
        """Shortcut for self.compute()."""
        return self.compute()
        
    def __str__(self):
        """
        String summary of the composite model.
        Lists each component and its immersion/emersion/opacity.
        """
        lines = ["CompositeModel with {} components:".format(len(self.components))]

        for name, model in self.components:
            p = model.params

            im = p.get("immersion") or p.get("immersion1")
            em = p.get("emersion")  or p.get("emersion1")
            op = p.get("opacity")        or p.get("opacity1")

            lines.append(f"  [{name}] {model.name}")
            if im is not None:
                lines.append(f"      immersion = {self.lightcurve.tref + im*u.s}")
            if em is not None:
                lines.append(f"      emersion  = {self.lightcurve.tref + em*u.s}")
            if op is not None:
                lines.append(f"      opacity        = {op}")

        return "\n".join(lines)


# Helper: attach models to LightCurve
def attach_to_lightcurve_class(LightCurveClass):
    """
    Attach model constructors as LightCurve methods.

    Examples
    --------
    >>> lc.square_well_model(immersion=..., emersion=..., opacity=0.8)
    >>> multi = lc.composite_model()
    """

    def _square_well_model(self, immersion, emersion, opacity=1.0
                           ):
        """Return a SquareWellModel linked to this LightCurve."""
        return SquareWellModel(
            lightcurve=self,
            immersion=immersion,
            emersion=emersion,
            opacity=opacity)
    def _double_square_well_model(self,
                                immersion1, emersion1, opacity1,
                                immersion2, emersion2, opacity2
                                ):
        """Return a DoubleSquareWellModel linked to this LightCurve."""
        return DoubleSquareWellModel(
            lightcurve=self,
            immersion1=immersion1, emersion1=emersion1, opacity1=opacity1,
            immersion2=immersion2, emersion2=emersion2, opacity2=opacity2)

    setattr(LightCurveClass, "double_square_well_model", _double_square_well_model)


    def _composite_model(self):
        """Return a new CompositeModel linked to this LightCurve."""
        model = CompositeModel()
        model.lightcurve = self
        return model

    setattr(LightCurveClass, "square_well_model", _square_well_model)
    setattr(LightCurveClass, "composite_model", _composite_model)


### Fast models
def occ_model(time, immersion, emersion, opacity,
              lambda_0, delta_lambda, distance, vel, exptime, d_star,
              npt_star=12, time_resolution_factor=10, flux_min=0.0, flux_max=1.0):
    model = SquareWellModel(
            lightcurve=None,
            immersion=immersion,
            emersion=emersion,
            opacity=opacity,
            npt_star=npt_star,
            lambda_0=lambda_0,
            delta_lambda=delta_lambda,
            distance=distance, 
            vel=vel, 
            exptime=exptime,
            d_star=d_star,
            time_resolution_factor=time_resolution_factor,
            flux_min=flux_min,
            flux_max=flux_max,
        )
    return model.compute(time)

def occ_model_double(time, immersion1, emersion1, opacity1,
                     immersion2, emersion2, opacity2,
                     lambda_0, delta_lambda, distance, vel, exptime, d_star,
                     npt_star=12, time_resolution_factor=10,
                     flux_min=0.0, flux_max=1.0):    
    model = DoubleSquareWellModel(
        lightcurve=None, 
        immersion1=immersion1,
        emersion1=emersion1,
        opacity1=opacity1,
        immersion2=immersion2,
        emersion2=emersion2,
        opacity2=opacity2,        
        npt_star=npt_star,
        lambda_0=lambda_0,
        delta_lambda=delta_lambda,
        distance=distance, 
        vel=vel, 
        exptime=exptime,
        d_star=d_star,
        time_resolution_factor=time_resolution_factor,
        flux_min=flux_min,
        flux_max=flux_max,
        )
    return model.compute(time=time)


    

def _lambda_weights(lambda_0, delta_lambda, n_lambda, response):
    """
    Optimized wavelength sampling for Fresnel integration.
        
    Assumes response is ALWAYS a tuple:
        response = (lamb_array, tlamb_array)
    """

    from scipy.interpolate import interp1d
    from numpy.polynomial.legendre import leggauss

    if response is None:
        if n_lambda == 1:
            return np.array([lambda_0]), np.array([1.0])
        lam_min = lambda_0 - delta_lambda/2
        lam_max = lambda_0 + delta_lambda/2
        lam_array = np.linspace(lam_min, lam_max, n_lambda)
        weights = np.ones_like(lam_array) / n_lambda
        weights /= weights.sum()
        return lam_array, weights

    lamb_resp = np.asarray(response[0])
    t_resp    = np.asarray(response[1])

    order = np.argsort(lamb_resp)
    lamb_resp = lamb_resp[order]
    t_resp    = t_resp[order]

    t_resp = np.clip(t_resp, 0, None)
    total = t_resp.sum()

    if total <= 0:
        raise ValueError("Spectral response has zero total throughput.")
    
    t_resp /= total

    if len(lamb_resp) <= 10:
        return lamb_resp, t_resp
    
    interp = interp1d(lamb_resp, t_resp, kind='linear',
                      bounds_error=False, fill_value=0.0)

    lam_min = lamb_resp.min()
    lam_max = lamb_resp.max()

    N = 5
    if n_lambda > 5:
        N = n_lambda
    
    xi, wi = leggauss(N)  

    lam_nodes = 0.5*(lam_max - lam_min)*xi + 0.5*(lam_max + lam_min)
    weights   = interp(lam_nodes) * wi

    total_w = weights.sum()
    if not np.isfinite(total_w) or total_w <= 0:
        raise RuntimeError("Invalid spectral weights derived from response.")

    weights /= total_w
    return lam_nodes, weights

def _plot_model(ax, lightcurve, flux_model,
                model_geometric=None,
                model_fresnel=None,
                model_star=None,
                time_model=None,
                show_components=False,
                title=None):

    ax = ax or plt.gca()

    if model_geometric is not None and time_model is not None:
        ax.plot(time_model, model_geometric, 'b-', label='Geometric', zorder=1)

        if show_components and time_model is not None:
            if model_fresnel is not None:
                ax.plot(time_model, model_fresnel, 'c-', label='Fresnel', zorder=1)
            if model_star is not None:
                ax.plot(time_model, model_star, 'g-', label='Star diam.', zorder=1)
    
    ax.plot(lightcurve.time, flux_model, 'r-', lw=1.0, label='Model', zorder=3)
    ax.scatter(lightcurve.time, flux_model, s=50, facecolors='none', edgecolors='r', zorder=3)

    if hasattr(lightcurve, "flux"):
        ax.plot(lightcurve.time, lightcurve.flux, "k.-", alpha=1, label='Observed')

    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Relative flux")
    ax.set_title(title or lightcurve.name)
    ax.legend()