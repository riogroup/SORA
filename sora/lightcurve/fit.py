"""
fit.py — Fitting routines for LightCurve
========================================

Preserves the philosophy and arguments of `LightCurve.occ_lcfit()` while exposing
a modern, concise interface:

    >>> model, chi2 = lc.fit(immersion_time=..., emersion_time=..., ...)

Methods supported:
    - 'chisqr'  : classic Monte Carlo (single-process)
    - 'fastchi' : Monte Carlo parallelizado (multiprocessing)
    - 'least_squares'/'ls' : Levenberg–Marquardt
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
from .model import occ_model

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

    def run(self, **kwargs):
        """
        Execute the fit preserving the occ_lcfit philosophy & args.
        Returns (model, ChiSquare).
        """
        allowed_kwargs = ['tmin', 'tmax', 'flux_min', 'flux_max',
                          'immersion_time', 'emersion_time', 'opacity',
                          'delta_t', 'dopacity', 'sigma', 'loop', 'verbose',
                          'sigma_result', 'method', 'threads', 'sigma_model']
        input_tests.check_kwargs(kwargs, allowed_kwargs=allowed_kwargs)

        if not hasattr(self.lc, 'flux'):
            raise ValueError('Fit curve is only possible when a LightCurve is instantiated with time and flux.')

        tmax = self.lc.time.max()
        tmin = self.lc.time.min()

        prelim = self.lc.occ_detect(tmin=tmin, tmax=tmax)

        delta_t = 2*self.lc.cycle
        loop = kwargs.get('loop', 10000)
        verbose = kwargs.get('verbose', True)

        immersion_time = tmin - self.lc.exptime
        do_immersion = False
        emersion_time = tmax + self.lc.exptime
        do_emersion = False
        opacity = kwargs.get('opacity', 1.0)
        delta_opacity = kwargs.get('dopacity', 0.0)
        do_opacity = ('dopacity' in kwargs)
        if do_opacity:
            warnings.warn("Fitting Opacity will be removed in future version")

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

        tmax = kwargs.get('tmax', tmax)
        tmin = kwargs.get('tmin', tmin)
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
                mdl = occ_model(time_m, t_i[i], t_e[i], opas[i],
                                self.lc.lambda_0, self.lc.delta_lambda,
                                self.lc.dist, self.lc.vel, self.lc.exptime, self.lc.d_star,
                                npt_star=12, time_resolution_factor=10,
                                flux_min=flux_min, flux_max=flux_max)
                chi2[i] = np.sum(((flux_m - mdl)**2) / (sig_m**2 + sigma_model**2))

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
            np.repeat(True, len(self.lc.flux)),
            flux_min=flux_min, flux_max=flux_max
        )

        self.lc.lc_sigma = sigma
        self.lc.chisquare = chisquare
        self.lc.opacity = opacity

        # Registro de resultados múltiplos no LightCurve
        # TODO [future feature]:
        # Integrar parâmetro `target` (Body, Ring, Satellite, etc.) aos métodos fit()
        # para vincular cada ajuste ao objeto físico correspondente.
        # Isso permitirá usar obj, obj.rings[...] e obj.satellites[...] como destino do fit

        if not hasattr(self.lc, "models"):
            self.lc.models = {}
        if not hasattr(self.lc, "chi2_maps"):
            self.lc.chi2_maps = {}
        if not hasattr(self.lc, "_fit_results"):
            self.lc._fit_results = {}

        label = f"fit{len(self.lc.models) + 1}"

        self.lc.models[label] = model
        self.lc.chi2_maps[label] = chisquare

        try:
            self.lc._fit_results[label] = {
                "type": "SquareWell",
                "immersion_time": self.lc._immersion,
                "immersion_err": getattr(self.lc, "immersion_err", np.nan),
                "emersion_time": self.lc._emersion,
                "emersion_err": getattr(self.lc, "emersion_err", np.nan),
                "opacity": opacity,
                "baseflux": flux_max,
                "bottomflux": flux_min,
                "curve_sigma": sigma.mean(),
            }
        except Exception as e:
            print(f"[WARN] Could not register fit result: {e}")
        return model, chisquare



# -----------------------------------------------------------
# Error function used by LS/DE (idêntico à filosofia do original)
# -----------------------------------------------------------
def _fit_error(parameters, time, flux, dflux, flux_min, flux_max,
               lambda_0, delta_lambda, distance, velocity, exptime, d_star,
               time_resolution_factor, npt_star):
    v = parameters.valuesdict()
    mdl = occ_model(time, v['immersion_time'], v['emersion_time'], v['opacity'],
                    lambda_0, delta_lambda, distance, velocity, exptime, d_star,
                    npt_star=npt_star, time_resolution_factor=time_resolution_factor,
                    flux_min=flux_min, flux_max=flux_max)
    return (flux - mdl)**2 / (dflux**2)


# -----------------------------------------------------------
# Public interface to be bound into LightCurve
# -----------------------------------------------------------
def fit(self, **kwargs):
    return _FitHandler(self).run(**kwargs)
