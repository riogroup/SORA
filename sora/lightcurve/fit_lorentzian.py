
"""
fit_lorentzian.py — Fitting routines for LorentzianModel
===========================================================

Matches LightCurve.fit() philosophy, but fits a Lorentzian dip:

    F(t) = flux_max - depth*(flux_max - flux_min)*P(t)

Parameters fitted:
    - center_time (s from tref)
    - fwhm (s)
    - depth (dimensionless; central drop in units of (flux_max - flux_min))

Supported methods:
    - 'chisqr'
    - 'fastchi'
    - 'least_squares'/'ls'
    - 'differential_evolution'/'de'

Returns:
    (model, ChiSquare)
"""

import numpy as np
import warnings
from multiprocessing import Pool

import astropy.units as u

from sora.extra import ChiSquare
from sora.stats import Parameters, least_squares, differential_evolution
from sora.config.visuals import progressbar
from sora.config import input_tests

from .model import LorentzianModel


def _mc_worker_lorentz(time, flux, sigma, bestchi,
                       center_time, fwhm, delta_t, delta_fwhm,
                       depth, delta_depth,
                       exptime, time_resolution_factor, flux_min, flux_max,
                       loop, verbose):

    rng = np.random.RandomState()

    c = center_time + delta_t * (2 * rng.random(loop) - 1)

    fw = fwhm + delta_fwhm * (2 * rng.random(loop) - 1)
    fw = np.clip(fw, 1e-6, np.inf)

    d = depth + delta_depth * (2 * rng.random(loop) - 1)
    d = np.clip(d, 0.0, 1.5)

    if bestchi is not None and loop > 0:
        c[0] = center_time
        fw[0] = fwhm
        d[0] = depth

    chi = np.empty(loop, dtype=float)

    iterator = progressbar(range(loop), "LightCurve Lorentz fit:") if verbose else range(loop)

    for i in iterator:
        mdl = LorentzianModel(
            lightcurve=None,
            center=c[i],
            fwhm=fw[i],
            depth=d[i],
            exptime=exptime,
            time_resolution_factor=time_resolution_factor,
            flux_min=flux_min,
            flux_max=flux_max,
        )
        m = mdl.compute(time=time)
        chi[i] = np.sum(((flux - m) ** 2) / (sigma ** 2))

    return [chi, c, fw, d]


def _fit_error_lorentz(parameters, time, flux, dflux,
                      exptime, time_resolution_factor, flux_min, flux_max):
    v = parameters.valuesdict()

    mdl = LorentzianModel(
        lightcurve=None,
        center=v["center_time"],
        fwhm=v["fwhm"],
        depth=v["depth"],
        exptime=exptime,
        time_resolution_factor=time_resolution_factor,
        flux_min=flux_min,
        flux_max=flux_max,
    )
    m = mdl.compute(time=time)
    return (flux - m) ** 2 / (dflux ** 2)


class _FitHandlerLorentz:
    def __init__(self, lc):
        self.lc = lc

    def run(self, clear_fits=True, **kwargs):
        allowed_kwargs = [
            "tmin", "tmax", "flux_min", "flux_max",
            "center_time", "fwhm", "depth",
            "delta_t", "delta_fwhm", "ddepth",
            "sigma", "loop", "verbose",
            "sigma_result", "method", "threads", "sigma_model",
            "time_resolution_factor",
        ]
        input_tests.check_kwargs(kwargs, allowed_kwargs=allowed_kwargs)

        if clear_fits and hasattr(self.lc, "clear_fits"):
            self.lc.clear_fits()

        if not hasattr(self.lc, "flux"):
            raise ValueError("Fit is only possible when LightCurve has time and flux.")

        method = str(kwargs.get("method") or "chisqr").lower()
        if method not in ["chisqr", "least_squares", "ls", "fastchi", "differential_evolution", "de"]:
            warnings.warn(f"Invalid method `{method}` provided. Setting to default.")
            method = "chisqr"

        tmax = kwargs.get("tmax", self.lc.time.max())
        tmin = kwargs.get("tmin", self.lc.time.min())
        mask = (self.lc.time >= tmin) & (self.lc.time <= tmax)

        prelim = self.lc.occ_detect(tmin=tmin, tmax=tmax)

        # Baseline / scale (same spirit as fit.py)
        flux_min = kwargs.get("flux_min", 1 - prelim["depth"])
        flux_max = kwargs.get("flux_max", prelim["baseline"])

        # Initial guesses from prelim detection
        imm0 = prelim["immersion_time"]
        eme0 = prelim["emersion_time"]
        dur0 = prelim["occultation_duration"]

        center_time = kwargs.get("center_time", 0.5 * (imm0 + eme0))
        fwhm = kwargs.get("fwhm", max(dur0, self.lc.exptime))
        # Depth as fraction of (flux_max - flux_min)
        # If prelim['depth'] is already a fractional drop (common in SORA), this is reasonable.
        depth = kwargs.get("depth", np.clip(prelim["depth"], 0.0, 1.0))

        # Ranges for MC sampling
        delta_t = kwargs.get("delta_t", 5 * prelim["time_err"])
        delta_fwhm = kwargs.get("delta_fwhm", 0.5 * fwhm)
        delta_depth = kwargs.get("ddepth", 0.1)

        loop = int(kwargs.get("loop", 10000))
        verbose = bool(kwargs.get("verbose", True))
        sigma_result = kwargs.get("sigma_result", 1.0)
        sigma_model = float(kwargs.get("sigma_model", 0.0))
        trf = float(kwargs.get("time_resolution_factor", 10.0))

        # Sigma handling (same behavior as fit.py)
        if "sigma" not in kwargs:
            if getattr(self.lc, "dflux", None) is not None:
                sigma = np.array(self.lc.dflux)
            else:
                sigma = "auto"
        else:
            s = kwargs["sigma"]
            if isinstance(s, (float, int)):
                sigma = np.repeat(float(s), len(self.lc.flux))
            elif s is None:
                sigma = np.array(self.lc.dflux) if getattr(self.lc, "dflux", None) is not None else "auto"
            else:
                sigma = np.array(s)

        if isinstance(sigma, str) and sigma == "auto":
            # estimate sigma outside the event window inferred by prelim
            mask_sigma = (((self.lc.time >= tmin) & (self.lc.time < imm0 - self.lc.exptime)) |
                          ((self.lc.time > eme0 + self.lc.exptime) & (self.lc.time <= tmax)))
            sig = self.lc.flux[mask_sigma].std(ddof=1)
            sigma = np.repeat(sig, len(self.lc.flux))

        time_m = self.lc.time[mask]
        flux_m = self.lc.flux[mask]
        sig_m = sigma[mask]

        # --- METHOD: chisqr (straight MC) ---
        if method == "chisqr":
            chi2 = np.full(loop, 999999.0, dtype=float)

            iterator = progressbar(range(loop), "LightCurve Lorentz fit:") if verbose else range(loop)

            c_s = center_time + delta_t * (2 * np.random.random(loop) - 1)
            fw_s = np.clip(fwhm + delta_fwhm * (2 * np.random.random(loop) - 1), 1e-6, np.inf)
            d_s = np.clip(depth + delta_depth * (2 * np.random.random(loop) - 1), 0.0, 1.5)

            for i in iterator:
                mdl = LorentzianModel(
                    lightcurve=None,
                    center=c_s[i],
                    fwhm=fw_s[i],
                    depth=d_s[i],
                    exptime=self.lc.exptime,
                    time_resolution_factor=trf,
                    flux_min=flux_min,
                    flux_max=flux_max,
                )
                m = mdl.compute(time=time_m)
                chi2[i] = np.sum(((flux_m - m) ** 2) / (sig_m ** 2 + sigma_model ** 2))

        # --- METHODS: ls / de / fastchi ---
        else:
            threads = int(kwargs.get("threads", 1))

            initial = Parameters()

            # center bounds
            initial.add(
                name="center_time",
                value=center_time,
                minval=center_time - delta_t,
                maxval=center_time + delta_t,
                free=True,
            )

            # fwhm bounds (positive)
            fmin = max(1e-6, fwhm - delta_fwhm)
            fmax = max(fmin * 1.001, fwhm + delta_fwhm)
            initial.add(
                name="fwhm",
                value=fwhm,
                minval=fmin,
                maxval=fmax,
                free=True,
            )

            # depth bounds
            dmin = max(0.0, depth - delta_depth)
            dmax = min(1.0, depth + delta_depth)
            initial.add(
                name="depth",
                value=depth,
                minval=dmin,
                maxval=dmax,
                free=True,
            )

            set_bestchi = False

            if method in ["least_squares", "ls"]:
                res = least_squares(
                    _fit_error_lorentz,
                    initial,
                    args=(time_m, flux_m, sig_m, self.lc.exptime, trf, flux_min, flux_max),
                    algorithm="trf",
                    sigma=sigma_result,
                )
                center_time = res.params["center_time"].value
                fwhm = res.params["fwhm"].value
                depth = res.params["depth"].value
                bestchi, set_bestchi = res.chisqr, True
                method = "fastchi"

            if method in ["differential_evolution", "de"]:
                res = differential_evolution(
                    _fit_error_lorentz,
                    initial,
                    args=(time_m, flux_m, sig_m, self.lc.exptime, trf, flux_min, flux_max),
                    sigma=sigma_result,
                )
                center_time = res.params["center_time"].value
                fwhm = res.params["fwhm"].value
                depth = res.params["depth"].value
                bestchi, set_bestchi = res.chisqr, True
                method = "fastchi"

            if method == "fastchi":
                bestchi = bestchi if set_bestchi else None

                per = int(np.ceil(loop / max(1, threads)))
                args_common = (
                    time_m, flux_m, sig_m, bestchi,
                    center_time, fwhm, delta_t, delta_fwhm,
                    depth, delta_depth,
                    self.lc.exptime, trf, flux_min, flux_max,
                )

                if verbose:
                    jobs = [(*args_common, per, True)] + [(*args_common, per, False) for _ in range(max(0, threads - 1))]
                else:
                    jobs = [(*args_common, per, False) for _ in range(threads)]

                with Pool(processes=threads) as pool:
                    results = pool.starmap(_mc_worker_lorentz, jobs)

                chi2, cc, ff, dd = [], [], [], []
                for r in results:
                    chi2.extend(r[0]); cc.extend(r[1]); ff.extend(r[2]); dd.extend(r[3])

                chi2 = np.array(chi2[:loop])
                cc = np.array(cc[:loop])
                ff = np.array(ff[:loop])
                dd = np.array(dd[:loop])

        # Build ChiSquare object
        chisquare = ChiSquare(
            chi2,
            len(self.lc.flux[mask]),
            center=cc if "cc" in locals() else None,
            fwhm=ff if "ff" in locals() else None,
            depth=dd if "dd" in locals() else None,
        )

        result_sigma = chisquare.get_nsigma(sigma=sigma_result)

        # Best-fit from nsigma
        if "center" in result_sigma:
            self.lc.center_time = result_sigma["center"][0]
            self.lc.center_err = result_sigma["center"][1]
        else:
            self.lc.center_err = None

        if "fwhm" in result_sigma:
            self.lc.fwhm = result_sigma["fwhm"][0]
            self.lc.fwhm_err = result_sigma["fwhm"][1]
        else:
            self.lc.fwhm_err = None

        if "depth" in result_sigma:
            self.lc.depth = result_sigma["depth"][0]
            self.lc.depth_err = result_sigma["depth"][1]
        else:
            self.lc.depth_err = None

        # Final model instance
        model = LorentzianModel(
            lightcurve=self.lc,
            center=center_time,
            fwhm=fwhm,
            depth=depth,
            time_resolution_factor=trf,
            flux_min=flux_min,
            flux_max=flux_max,
        )

        model.center_err = self.lc.center_err
        model.fwhm_err = self.lc.fwhm_err
        model.depth_err = self.lc.depth_err

        # Store results
        if not hasattr(self.lc, "models"):
            self.lc.models = {}
        if not hasattr(self.lc, "chi2_maps"):
            self.lc.chi2_maps = {}
        if not hasattr(self.lc, "_fit_results"):
            self.lc._fit_results = {}

        label = f"fit{len(self.lc.models) + 1}"

        self.lc.models[label] = model
        self.lc.chi2_maps[label] = chisquare

        center_time_abs = self.lc.tref + center_time * u.s if hasattr(self.lc, "tref") else center_time

        self.lc._fit_results[label] = {
            "type": "Lorentzian",
            "center_time": center_time_abs,
            "center_err": self.lc.center_err,
            "fwhm": self.lc.fwhm,
            "fwhm_err": self.lc.fwhm_err,
            "depth": self.lc.depth,
            "depth_err": self.lc.depth_err,
            "baseflux": flux_max,
            "bottomflux": flux_min,
            "curve_sigma": float(np.mean(sigma[mask])),
        }

        return model, chisquare


def fit_lorentzian(self, clear_fits=True, **kwargs):
    """
    Fit the LightCurve using a Lorentzian dip model.

    Returns (model, ChiSquare).
    """
    return _FitHandlerLorentz(self).run(clear_fits=clear_fits, **kwargs)
