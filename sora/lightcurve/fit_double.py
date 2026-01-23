"""
fit_double.py — Fitting routines for LightCurve (Double Square Well)
====================================================================

Extends LightCurve.occ_lcfit() to handle two consecutive square-well regions coherently.

Usage example
-------------
    >>> model, chi2 = lc.fit_double(
    ...     immersion1=..., emersion1=..., opacity1=...,
    ...     immersion2=..., emersion2=..., opacity2=...,
    ...     method='fastchi', loop=5000, threads=4, verbose=True
    ... )

Methods supported:
    - 'chisqr'  : classic Monte Carlo (single-process)
    - 'fastchi' : Monte Carlo parallelizado (multiprocessing)
    - 'least_squares'/'ls' : Levenberg–Marquardt
    - 'differential_evolution'/'de' : DE global

Returns:
    (model, ChiSquare)
"""

import numpy as np
from multiprocessing import Pool
import astropy.units as u

from sora.extra import ChiSquare
from sora.stats import Parameters, least_squares, differential_evolution
from sora.config.visuals import progressbar
from sora.config import input_tests
from .model import *


# Internal worker for parallel fastchi
def _mc_worker(time, flux, sigma, bestchi,
               immersion1, emersion1, opacity1,
               immersion2, emersion2, opacity2,
               delta_t, delta_opacity,
               lambda_0, delta_lambda, distance, vel, exptime, d_star,
               npt_star, time_resolution_factor, flux_min, flux_max,
               loop, verbose):

    rng = np.random.RandomState()

    im1 = immersion1 + delta_t * (2*rng.random(loop) - 1)
    em1 = emersion1  + delta_t * (2*rng.random(loop) - 1)
    im2 = immersion2 + delta_t * (2*rng.random(loop) - 1)
    em2 = emersion2  + delta_t * (2*rng.random(loop) - 1)
    op1 = np.clip(opacity1 + delta_opacity * (2*rng.random(loop) - 1), 0.0, 1.0)
    op2 = np.clip(opacity2 + delta_opacity * (2*rng.random(loop) - 1), 0.0, 1.0)

    if bestchi is not None and loop > 0:
        im1[0], em1[0], op1[0] = immersion1, emersion1, opacity1
        im2[0], em2[0], op2[0] = immersion2, emersion2, opacity2

    chi = np.empty(loop, dtype=float)

    if verbose:
        iterator = progressbar(range(loop), 'LightCurve double fit:')
    else:
        iterator = range(loop)

    swap = im1 > em1
    im1[swap], em1[swap] = em1[swap], im1[swap]
    swap = im2 > em2
    im2[swap], em2[swap] = em2[swap], im2[swap]

    for i in iterator:
        mdl = DoubleSquareWellModel(
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
        chi[i] = np.sum(((flux - mdl.compute(time=time))**2) / (sigma**2))

    return [chi, im1, em1, op1, im2, em2, op2]


# Error function for LS/DE
def _fit_error_double(parameters, time, flux, dflux,
                      flux_min, flux_max,
                      lambda_0, delta_lambda, distance, vel,
                      exptime, d_star,
                      time_resolution_factor, npt_star):
    v = parameters.valuesdict()
    mdl = DoubleSquareWellModel(
        lightcurve=None, 
        immersion1=v['immersion1'],
        emersion1=v['emersion1'],
        opacity1=v['opacity1'],
        immersion2=v['immersion2'],
        emersion2=v['emersion2'],
        opacity2=v['opacity2'],        
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
    
    return (flux - mdl.compute(time))**2 / (dflux**2)


# LightCurve handler
class _FitDoubleHandler:
    def __init__(self, lc):
        self.lc = lc

    def run(self, **kwargs):
        allowed_kwargs = ['tmin', 'tmax', 'flux_min', 'flux_max',
                          'immersion1', 'emersion1', 'opacity1',
                          'immersion2', 'emersion2', 'opacity2',
                          'delta_t', 'dopacity', 'sigma', 'loop', 'verbose',
                          'sigma_result', 'method', 'threads', 'sigma_model', 
                          'feature']
        input_tests.check_kwargs(kwargs, allowed_kwargs=allowed_kwargs)

        if not hasattr(self.lc, 'flux'):
            raise ValueError("LightCurve must have 'time' and 'flux' for fitting.")

        tmin = kwargs.get('tmin', self.lc.time.min())
        tmax = kwargs.get('tmax', self.lc.time.max())

        required = ['immersion1','emersion1','opacity1','immersion2','emersion2','opacity2']
        missing = [r for r in required if r not in kwargs]
        if missing:
            raise ValueError(f"Missing initial parameters for double fit: {', '.join(missing)}")

        im1, em1, op1 = kwargs['immersion1'], kwargs['emersion1'], kwargs['opacity1']
        im2, em2, op2 = kwargs['immersion2'], kwargs['emersion2'], kwargs['opacity2']

        sigma_model = kwargs.get('sigma_model', 0.0)
        delta_t = kwargs.get('delta_t', 2*self.lc.cycle)
        delta_opacity = kwargs.get('dopacity', 0.0)
        loop = int(kwargs.get('loop', 10000))

        t_i1 = im1 + delta_t*(2*np.random.random(loop)-1)
        t_e2 = em1 + delta_t*(2*np.random.random(loop)-1)
        op_1 = np.clip(op1 + delta_opacity*(2*np.random.random(loop)-1), 0, 1)
        t_i2 = im2 + delta_t*(2*np.random.random(loop)-1)
        t_e2 = em2 + delta_t*(2*np.random.random(loop)-1)
        op_2 = np.clip(op2 + delta_opacity*(2*np.random.random(loop)-1), 0, 1)

        feature = str(kwargs.get('feature', 'body')).lower().strip()
        if feature in ('main', 'primary', 'object', 'limb'):
            feature = 'body'
        if feature in ('ring', 'rings'):
            feature = 'ring'
        
        verbose = bool(kwargs.get('verbose', True))
        sigma_result = kwargs.get('sigma_result', 1)
        threads = int(kwargs.get('threads', 1))
        method = str(kwargs.get('method', 'chisqr')).lower()
        if method not in ['chisqr', 'least_squares', 'ls', 'fastchi', 'differential_evolution', 'de']:
            warnings.warn(f'Invalid method `{method}` provided. Setting to default.')
            method = 'chisqr'

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
            mask_sigma = (((self.lc.time >= tmin) & (self.lc.time < im1 - self.lc.exptime)) |
                          ((self.lc.time > em2 + self.lc.exptime) & (self.lc.time <= tmax)))
            sig = self.lc.flux[mask_sigma].std(ddof=1)
            sigma = np.repeat(sig, len(self.lc.flux))

        flux_min = kwargs.get('flux_min', 0)
        flux_max = kwargs.get('flux_max', 1)

        bestchi, set_bestchi = None, False

        if method in ['least_squares', 'ls', 'differential_evolution', 'de']:
            initial = Parameters()

            for tag, val in zip(
                ['immersion1','emersion1','immersion2','emersion2'],
                [im1, em1, im2, em2]
            ):
                initial.add(name=tag, value=val,
                            minval=val - delta_t, maxval=val + delta_t, free=True)

            for tag, val in zip(['opacity1','opacity2'], [op1, op2]):
                initial.add(name=tag, value=val,
                            minval=max(0, val - delta_opacity),
                            maxval=min(1, val + delta_opacity),
                            free=(delta_opacity > 0))

            args = (self.lc.time[mask], self.lc.flux[mask], 
                    sigma[mask], flux_min, flux_max,
                    self.lc.lambda_0, self.lc.delta_lambda, 
                    self.lc.dist, self.lc.vel,
                    self.lc.exptime, self.lc.d_star, 10, 12)

            if method in ['least_squares', 'ls']:
                res = least_squares(_fit_error_double, initial, args=args,
                                    algorithm='trf', sigma=sigma_result)
            else:
                res = differential_evolution(_fit_error_double, initial, args=args,
                                             sigma=sigma_result)

            im1 = res.params['immersion1'].value
            em1 = res.params['emersion1'].value
            op1 = res.params['opacity1'].value
            im2 = res.params['immersion2'].value
            em2 = res.params['emersion2'].value
            op2 = res.params['opacity2'].value
            bestchi, set_bestchi = res.chisqr, True
            method = 'fastchi'  # refinement

        if method in ['chisqr', 'fastchi']:
            if method == 'chisqr':
                chi2 = np.empty(loop, float)
                rng = progressbar(range(loop), 'LightCurve double fit:') if verbose else range(loop)
                for i in rng:
                    mdl = DoubleSquareWellModel(
                        self.lc.time[mask],
                        t_i1, t_e2, op_1,
                        t_i2, t_e2, op_2,                                           
                        self.lc.lambda_0, self.lc.delta_lambda,
                        self.lc.dist, self.lc.vel, 
                        self.lc.exptime, self.lc.d_star,
                        npt_star=12, time_resolution_factor=10,
                        flux_min=flux_min, flux_max=flux_max
                        )
                    chi2[i] = np.sum(((lc.flux[mask] - mdl.compute(timeself.lc.time[mask]))**2) / (sigma[mask]**2 + sigma_model**2))

                im1s = np.repeat(im1, loop); em1s = np.repeat(em1, loop)
                op1s = np.repeat(op1, loop); im2s = np.repeat(im2, loop)
                em2s = np.repeat(em2, loop); op2s = np.repeat(op2, loop)

            else: 
                per = int(np.ceil(loop / max(1, threads)))
                args_common = (self.lc.time[mask], self.lc.flux[mask], sigma[mask], bestchi,
                               im1, em1, op1, im2, em2, op2,
                               delta_t, delta_opacity,
                               self.lc.lambda_0, self.lc.delta_lambda, 
                               self.lc.dist, self.lc.vel, self.lc.exptime, self.lc.d_star,
                               12, 10, flux_min, flux_max)
                jobs = [(*args_common, per, verbose if i == 0 else False) for i in range(threads)]
                with Pool(processes=threads) as pool:
                    results = pool.starmap(_mc_worker, jobs)

                chi2, im1s, em1s, op1s, im2s, em2s, op2s = [], [], [], [], [], [], []
                for r in results:
                    chi2.extend(r[0]); im1s.extend(r[1]); em1s.extend(r[2]); op1s.extend(r[3])
                    im2s.extend(r[4]); em2s.extend(r[5]); op2s.extend(r[6])

                chi2 = np.array(chi2[:loop])
                im1s, em1s, op1s = np.array(im1s[:loop]), np.array(em1s[:loop]), np.array(op1s[:loop])
                im2s, em2s, op2s = np.array(im2s[:loop]), np.array(em2s[:loop]), np.array(op2s[:loop])

        chisquare = ChiSquare(chi2, len(self.lc.flux[mask]),
                              immersion1=im1s, emersion1=em1s, opacity1=op1s,
                              immersion2=im2s, emersion2=em2s, opacity2=op2s)

        result_sigma = chisquare.get_nsigma(sigma=sigma_result)
        for key in ['immersion1','emersion1','opacity1','immersion2','emersion2','opacity2']:
            if key in result_sigma:
                locals()[key[:2] + (key[2:] if len(key)>2 else '')] = result_sigma[key][0]

        model = DoubleSquareWellModel(
            lightcurve=self.lc,
            immersion1=im1, emersion1=em1, opacity1=op1,
            immersion2=im2, emersion2=em2, opacity2=op2
        )
        model.compute(time=self.lc.time)

        self.lc.chisquare = chisquare
        self.lc.double_best_params = dict(
            immersion1=im1, emersion1=em1, opacity1=op1,
            immersion2=im2, emersion2=em2, opacity2=op2
        )

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
                "type": "DoubleSquareWell",
                "immersion1_time": self.lc.tref + im1 * u.s,
                "emersion1_time": self.lc.tref + em1 * u.s,
                "opacity1": op1,
                "immersion2_time": self.lc.tref + im2 * u.s,
                "emersion2_time": self.lc.tref + em2 * u.s,
                "opacity2": op2,
                "baseflux": flux_max,
                "bottomflux": flux_min,
                "curve_sigma": sigma.mean(),
                "feature": feature
                }
        except Exception:
            pass

        return model, chisquare

def fit_double(self, **kwargs):
    return _FitDoubleHandler(self).run(**kwargs)
