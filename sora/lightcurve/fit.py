"""
fit.py — Fitting routines for LightCurve
========================================

Preserves the philosophy and arguments of `LightCurve.occ_lcfit()` while exposing
a modern, concise interface:

    >>> model, chi2 = lc.fit(immersion_time=..., emersion_time=..., ...)

Methods supported:
    - 'chisqr'  : classic Monte Carlo (single-process)
    - 'fastchi' : Monte Carlo parallelizado (multiprocessing)
    - 'least_squares'/'ls' : Levenberg-Marquardt
    - 'differential_evolution'/'de' : DE global

Returns:
    (model, ChiSquare)
"""

import numpy as np
import warnings
from multiprocessing import Pool

from sora.extra import ChiSquare
from sora.stats import Parameters, least_squares, differential_evolution
from sora.config.visuals import progressbar
from sora.config import input_tests
from .model import *

import astropy.units as u  

# Internal worker for parallel fastchi
def _mc_worker(time, flux, sigma, bestchi,
               immersion_time, emersion_time, delta_t,
               opacity, delta_opacity,
               lambda_0, delta_lambda, dist, vel, exptime, d_star,
               npt_star, time_resolution_factor, flux_min, flux_max,
               loop, verbose):

    rng = np.random.RandomState()

    im = immersion_time + delta_t * (2*rng.random(loop) - 1)
    em = emersion_time + delta_t * (2*rng.random(loop) - 1)
    op = opacity + delta_opacity * (2*rng.random(loop) - 1)
    op = np.clip(op, 0.0, 1.0)

    if bestchi is not None and loop > 0:
        im[0] = immersion_time
        em[0] = emersion_time
        op[0] = opacity

    chi = np.empty(loop, dtype=float)

    if verbose:
        iterator = progressbar(range(loop), 'LightCurve fit:')
    else:
        iterator = range(loop)

    for i in iterator:
        mdl = occ_model(time, im[i], em[i], op[i],
                        lambda_0, delta_lambda, dist, vel, exptime, d_star,
                        npt_star=npt_star, time_resolution_factor=time_resolution_factor,
                        flux_min=flux_min, flux_max=flux_max)
        chi[i] = np.sum(( (flux - mdl)**2 ) / (sigma**2))

    return [chi, im, em, op]


# LightCurve-facing handler
class _FitHandler:
    def __init__(self, lc):
        self.lc = lc

    def run(self, clear_fits=True, **kwargs):
        """
        Execute the fit preserving the occ_lcfit philosophy & args.
        Returns (model, ChiSquare).
        """
        allowed_kwargs = ['tmin', 'tmax', 'flux_min', 'flux_max',
                          'immersion_time', 'emersion_time', 'opacity',
                          'delta_t', 'dopacity', 'sigma', 'loop', 'verbose',
                          'sigma_result', 'method', 'threads', 'sigma_model']
        input_tests.check_kwargs(kwargs, allowed_kwargs=allowed_kwargs)

        if clear_fits and hasattr(self.lc, "clear_fits"):
            self.lc.clear_fits()

        if not hasattr(self.lc, 'flux'):
            raise ValueError('Fit curve is only possible when a LightCurve is instantiated with time and flux.')

        tmax = kwargs.get('tmax', self.lc.time.max())
        tmin = kwargs.get('tmin', self.lc.time.min())

        prelim = self.lc.occ_detect(tmin=tmin, tmax=tmax)

        delta_t = 2*self.lc.cycle
        loop = kwargs.get('loop', 10000)
        verbose = kwargs.get('verbose', True)

        immersion_time = self.lc.time.min() - self.lc.exptime
        do_immersion = False
        emersion_time = self.lc.time.max() + self.lc.exptime
        do_emersion = False
        opacity = kwargs.get('opacity', 1.0)
        delta_opacity = kwargs.get('dopacity', 0.0)
        do_opacity = ('dopacity' in kwargs)

        if ('immersion_time' not in kwargs) and ('emersion_time' not in kwargs):
            immersion_time = prelim['immersion_time']
            do_immersion = True
            emersion_time = prelim['emersion_time']
            do_emersion = True
            delta_t = 5*prelim['time_err']
            tmax = emersion_time + 2*prelim['occultation_duration']
            tmin = immersion_time - 2*prelim['occultation_duration']
            if 2*prelim['occultation_duration'] < 10*self.lc.cycle:
                tmax = emersion_time + 10*self.lc.cycle
                tmin = immersion_time - 10*self.lc.cycle


        sigma_model = kwargs.get('sigma_model', 0.0)
        delta_t = kwargs.get('delta_t', delta_t)

        if 'immersion_time' in kwargs:
            immersion_time = kwargs['immersion_time']
            do_immersion = True
        if 'emersion_time' in kwargs:
            emersion_time = kwargs['emersion_time']
            do_emersion = True

        t_i = immersion_time + delta_t*(2*np.random.random(loop) - 1)
        t_e = emersion_time + delta_t*(2*np.random.random(loop) - 1)

        mask = (self.lc.time >= tmin) & (self.lc.time <= tmax)

        if 'sigma' not in kwargs:
            if getattr(self.lc, 'dflux', None) is not None:
                sigma = np.array(self.lc.dflux)
            else:
                sigma = 'auto'
        else:
            if isinstance(kwargs['sigma'], (float, int)):
                sigma = np.repeat(kwargs['sigma'], len(self.lc.flux))
            elif kwargs['sigma'] is None:
                sigma = np.array(self.lc.dflux) if getattr(self.lc, 'dflux', None) is not None else 'auto'
            else:
                sigma = np.array(kwargs['sigma'])

        if isinstance(sigma, str) and (sigma == 'auto'):
            mask_sigma = (((self.lc.time >= tmin) & (self.lc.time < immersion_time - self.lc.exptime)) |
                          ((self.lc.time > emersion_time + self.lc.exptime) & (self.lc.time <= tmax)))
            sig = self.lc.flux[mask_sigma].std(ddof=1)
            sigma = np.repeat(sig, len(self.lc.flux))

        # opacidades para MC
        opas = opacity + delta_opacity*(2*np.random.random(loop) - 1)
        opas = np.clip(opas, 0.0, 1.0)

        # baseline e bottom (iguais à heurística do occ_lcfit)
        flux_min = kwargs.get('flux_min', 1 - prelim['depth'])
        flux_max = kwargs.get('flux_max', prelim['baseline'])
        sigma_result = kwargs.get('sigma_result', 1)

        # garantir ordem (t_i <= t_e)
        swap = t_i > t_e
        t_i[swap], t_e[swap] = t_e[swap], t_i[swap]

        # método
        method = str(kwargs.get('method') or 'chisqr').lower()
        if method not in ['chisqr', 'least_squares', 'ls', 'fastchi', 'differential_evolution', 'de']:
            warnings.warn(f'Invalid method `{method}` provided. Setting to default.')
            method = 'chisqr'

        # Dispatch
        if method == 'chisqr':
            chi2 = 999999*np.ones(loop)
            if verbose:
                rng = progressbar(range(loop), 'LightCurve fit:')
            else:
                rng = range(loop)

            time_m = self.lc.time[mask]
            flux_m = self.lc.flux[mask]
            sig_m = sigma[mask]

            for i in rng:
                mdl = SquareWellModel(
                                    t_i[i], t_e[i], opas[i],
                                    self.lc.lambda_0, self.lc.delta_lambda,
                                    self.lc.dist, self.lc.vel, self.lc.exptime, self.lc.d_star,
                                    npt_star=12, time_resolution_factor=10,
                                    flux_min=flux_min, flux_max=flux_max)
                chi2[i] = np.sum(((flux_m - mdl.compute(time=time_m))**2) / (sig_m**2 + sigma_model**2))

        else:
            threads = int(kwargs.get('threads', 1))

            initial = Parameters()

            im_vary = not ((not do_immersion) or (delta_t == 0))
            initial.add(name='immersion_time', value=(immersion_time if do_immersion else tmin),
                        minval=(-np.inf if not im_vary else (immersion_time - delta_t)),
                        maxval=( np.inf if not im_vary else (immersion_time + delta_t)),
                        free=im_vary)

            em_vary = not ((not do_emersion) or (delta_t == 0))
            initial.add(name='emersion_time', value=(emersion_time if do_emersion else tmax),
                        minval=(-np.inf if not em_vary else (emersion_time - delta_t)),
                        maxval=( np.inf if not em_vary else (emersion_time + delta_t)),
                        free=em_vary)

            opacity_vary = (do_opacity and (delta_opacity != 0))
            minop = 0 if (opacity - delta_opacity) < 0 else (opacity - delta_opacity)
            maxop = 1 if (opacity + delta_opacity) > 1 else (opacity + delta_opacity)
            initial.add(name='opacity', value=(opacity if do_opacity else 1.0),
                        minval=(-np.inf if not opacity_vary else minop),
                        maxval=( np.inf if not opacity_vary else maxop),
                        free=opacity_vary)

            if (not im_vary) and (not em_vary) and (not opacity_vary):
                raise SystemExit('No parameters are allowed to vary, please check your `LightCurve.fit` input.')

            set_bestchi = False

            # LS
            if method in ['least_squares', 'ls']:
                res = least_squares(_fit_error, initial,
                                    args=(self.lc.time[mask], self.lc.flux[mask], sigma[mask],
                                          flux_min, flux_max,
                                          self.lc.lambda_0, self.lc.delta_lambda,
                                          self.lc.dist, self.lc.vel, self.lc.exptime, self.lc.d_star,
                                          10, 12),
                                    algorithm='trf', sigma=sigma_result)

                immersion_time = res.params['immersion_time'].value
                emersion_time  = res.params['emersion_time'].value
                opacity        = res.params['opacity'].value
                bestchi, set_bestchi = res.chisqr, True
                method = 'fastchi'   # refina via MC paralelo

            # DE
            if method in ['differential_evolution', 'de']:
                res = differential_evolution(_fit_error, initial,
                                             args=(self.lc.time[mask], self.lc.flux[mask], sigma[mask],
                                                   flux_min, flux_max,
                                                   self.lc.lambda_0, self.lc.delta_lambda,
                                                   self.lc.dist, self.lc.vel, self.lc.exptime, self.lc.d_star,
                                                   10, 12),
                                             sigma=sigma_result)

                immersion_time = res.params['immersion_time'].value
                emersion_time  = res.params['emersion_time'].value
                opacity        = res.params['opacity'].value
                bestchi, set_bestchi = res.chisqr, True
                method = 'fastchi'   # refina via MC paralelo

            # FASTCHI
            if method == 'fastchi':
                bestchi = bestchi if set_bestchi else None
                time_m = self.lc.time[mask]
                flux_m = self.lc.flux[mask]
                sig_m  = sigma[mask]

                # distribuir carga por processos
                per = int(np.ceil(loop/max(1, threads)))
                args_common = (time_m, flux_m, sig_m, bestchi,
                               immersion_time, emersion_time, delta_t,
                               opacity, delta_opacity,
                               self.lc.lambda_0, self.lc.delta_lambda, self.lc.dist, self.lc.vel,
                               self.lc.exptime, self.lc.d_star, 12, 10, flux_min, flux_max)

                if verbose:
                    # 1 processo verboso + (threads-1) silenciosos
                    jobs = [( *args_common, per, True )]
                    for _ in range(max(0, threads-1)):
                        jobs.append( (*args_common, per, False) )
                else:
                    jobs = [ (*args_common, per, False) for _ in range(threads) ]

                with Pool(processes=threads) as pool:
                    results = pool.starmap(_mc_worker, jobs)

                # concatenar
                chi2, ti, te, oo = [], [], [], []
                for r in results:
                    chi2.extend(r[0]); ti.extend(r[1]); te.extend(r[2]); oo.extend(r[3])

                chi2 = np.array(chi2[:loop])
                t_i  = np.array(ti [:loop])
                t_e  = np.array(te [:loop])
                opas = np.array(oo [:loop])

        # construir ChiSquare
        kkwargs = {}
        if (do_immersion) and (delta_t > 0):
            kkwargs['immersion'] = t_i
        if (do_emersion) and (delta_t > 0):
            kkwargs['emersion'] = t_e
        if (do_opacity) and (delta_opacity > 0):
            kkwargs['opacity'] = opas

        chisquare = ChiSquare(chi2, len(self.lc.flux[mask]), **kkwargs)

        # melhores valores (como no occ_lcfit)
        result_sigma = chisquare.get_nsigma(sigma=sigma_result)

        if 'immersion' in result_sigma:
            self.lc._immersion = self.lc.tref + result_sigma['immersion'][0]*u.s
            self.lc.immersion_err = result_sigma['immersion'][1]
            immersion_time = result_sigma['immersion'][0]
        else:
            try:
                immersion_time = (self.lc._immersion.jd - self.lc.tref.jd)*u.d.to('s')
            except Exception:
                pass

        if 'emersion' in result_sigma:
            self.lc._emersion = self.lc.tref + result_sigma['emersion'][0]*u.s
            self.lc.emersion_err = result_sigma['emersion'][1]
            emersion_time = result_sigma['emersion'][0]
        else:
            try:
                emersion_time = (self.lc._emersion.jd - self.lc.tref.jd)*u.d.to('s')
            except Exception:
                pass

        if 'opacity' in result_sigma:
            opacity = result_sigma['opacity'][0]

        # roda o modelo final e salva no objeto
        model = self.lc.occ_model(
            immersion_time, emersion_time, opacity,
            flux_min=flux_min, flux_max=flux_max
        )

        self.lc.lc_sigma = sigma
        self.lc.chisquare = chisquare
        self.lc.opacity = opacity

        # Registro de resultados múltiplos no LightCurve
        # TODO [future feature]:
        # Integrar parâmetro `target` (Body, Ring, Satellite, etc.) aos métodos fit()
        # para vincular cada ajuste ao objeto físico correspondente.

        if not hasattr(self.lc, "models"):
            self.lc.models = {}
        if not hasattr(self.lc, "chi2_maps"):
            self.lc.chi2_maps = {}
        if not hasattr(self.lc, "_fit_results"):
            self.lc._fit_results = {}

        label = f"fit{len(self.lc.models) + 1}"

        nsig1 = chisquare.get_nsigma(sigma=sigma_result)
        model.immersion_err = nsig1["immersion"][1] if 'immersion' in nsig1.keys() else None
        model.emersion_err = nsig1["emersion"][1] if 'emersion' in nsig1.keys() else None
        model.opacity_err = nsig1["opacity"][1] if 'opacity' in nsig1.keys() else None

        self.lc.models[label] = model
        self.lc.chi2_maps[label] = chisquare

        self.lc._fit_results[label] = {
            "type": "SquareWell",
            "immersion_time": model.immersion,
            "immersion_err": model.immersion_err,
            "emersion_time": model.emersion,
            "emersion_err": model.emersion_err,
            "opacity": model.opacity,
            "opacity_err": model.opacity_err,
            "baseflux": flux_max,
            "bottomflux": flux_min,
            "curve_sigma": sigma.mean(),
        }

        return model, chisquare



# -----------------------------------------------------------
# Error function used by LS/DE (idêntico à filosofia do original)
# -----------------------------------------------------------
def _fit_error(parameters, time, flux, dflux, flux_min, flux_max,
               lambda_0, delta_lambda, distance, vel, exptime, d_star,
               time_resolution_factor, npt_star):
    v = parameters.valuesdict()
    model = SquareWellModel(
            lightcurve=None,
            immersion=v['immersion_time'],
            emersion= v['emersion_time'],
            opacity=v['opacity'],
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
    mdl=model.compute(time)
    return (flux - mdl)**2 / (dflux**2)

def fit(self, clear_fits=True, **kwargs):
    """ Fit the LightCurve using a square-well occultation model.
    The method returns a tuple:

        (model, ChiSquare)

    where `model` is a `SquareWellModel` instance evaluated with the best-fit
    parameters, and `ChiSquare` contains the sampled chi-square distribution
    and uncertainties.

    Parameters
    ----------
    clear_fits : bool, optional, default=True
        If True, all stored models and chi-square maps attached to the
        LightCurve are cleared before performing the new fit. Set to False
        to accumulate multiple labelled fits (`fit1`, `fit2`, ...).

    tmin : int, float, optional
        Minimum time (seconds relative to `tref`) to consider in the fit.
        Defaults to the earliest timestamp in the LightCurve.

    tmax : int, float, optional
        Maximum time (seconds relative to `tref`) to consider in the fit.
        Defaults to the latest timestamp in the LightCurve.

    flux_min : int, float, optional
        Bottom flux level used by the square-well model.
        If not provided, the routine estimates it from the preliminary
        occultation detection.

    flux_max : int, float, optional
        Baseline flux level. If not provided, it is estimated from the
        preliminary detection.

    immersion_time : int, float, optional
        Initial guess for immersion time (seconds from `tref`).
        If not given, a preliminary detection (`occ_detect`) is used.

    emersion_time : int, float, optional
        Initial guess for emersion time (seconds from `tref`).
        If not given, a preliminary detection is used.

    opacity : int, float, default=1.0
        Initial opacity. 1 = opaque, 0 = transparent.
        The opacity implicitly includes Fresnel diffraction effects and
        the cross-sectional attenuation expected from ring particles.

    delta_t : int, float, optional
        Half-range used to vary immersion and emersion during the Monte Carlo
        sampling. If not provided, defaults to twice the light-curve cycle or
        the preliminary-detection time uncertainty.

    dopacity : int, float, optional
        Half-range used to vary opacity during Monte Carlo sampling.

    sigma : int, float, array_like, or 'auto', optional
        Flux uncertainties. If None, uses `self.dflux`.
        If 'auto', sigma is estimated from regions outside the occultation.

    sigma_model : int, float, optional, default=0
        Additional model uncertainty (in flux units) added in quadrature to
        `sigma`. Used only with `chisqr`.

    loop : int, optional, default=10000
        Number of Monte Carlo trials.

    verbose : bool, optional, default=True
        If True, displays progress bars and diagnostics during sampling.

    sigma_result : int, float, optional
        Sigma level used to compute final parameter uncertainties from the
        ChiSquare object (e.g. sigma_result=1 for 1sigma errors).

    method : str, optional, default='chisqr'
        Fitting method. Accepted values are:
        - `chisqr`  : classic Monte Carlo chi-square sampling (single process)
        - `fastchi` : parallel Monte Carlo accelerated with multiprocessing
        - `least_squares` / `ls` :
              Levenberg-Marquardt minimization followed by a MC refinement
        - `differential_evolution` / `de` :
              global genetic optimization followed by a MC refinement

    threads : int, optional, default=1
        Number of worker processes used by the `fastchi` or DE refinement.
        Ignored for `chisqr`.

    Returns
    -------
    model : `SquareWellModel`
        Square-well model computed with the best-fit immersion, emersion
        and opacity.

    chi2 : `sora.extra.ChiSquare`
        ChiSquare object containing the Monte Carlo distribution, best-fit
        values, and parameter uncertainties.

    Notes
    -----
    - If neither `immersion_time` nor `emersion_time` is given, a preliminary
      detection is performed to estimate their initial values.

    - The fitting region (`tmin`, `tmax`) is automatically expanded around the
      event if the preliminary occultation duration is short.

    - All results are stored inside the LightCurve under labelled entries
      (`lc.models['fit1']`, `lc.chi2_maps['fit1']`, etc.).

    - The final model and uncertainties follow the same conventions as the
      original SORA `occ_lcfit()` routine, ensuring backward compatibility.
    """

    return _FitHandler(self).run(clear_fits=clear_fits, **kwargs)
