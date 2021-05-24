import warnings
import numpy as np

__all__ = ['occ_detect']


def occ_detect(flux, dflux, time, cycle, maximum_duration=None, dur_step=None, snr_limit=None,
               n_detections=None, plot=False):
    """Detects automatically the occultation event in the light curve.

    Detects a 'square well' shaped transit. All parameters are optional.


    Parameters
    ----------
    flux : `float` array
        Flux of the time series. Dependent variable.

    dflux: `float` array
        Error in the flux. Error in the dependent variable.

    time: `float` array
        Time variable. Independent variable.

    cycle: `float`
        Sampling value of the time series.

    maximum_duration : `float`, default: light curve time span
        Maximum duration of the occultation event.

    dur_step : `float`, default: 1/2 cycle
        Step size to sweep occultation duration event.

    snr_limit : `float`, default=None
        Minimum occultation SNR.

    n_detections : `int`, default=1
        Number of detections regardless the SNR. `n_detections` is superseded by
        `snr_limit`.

    plot : `boolean`, default=False
        True if output plots are desired.


    Returns
    -------
    OrderedDict : `dict`
        An ordered dictionary of :attr:`name`::attr:`value` pairs for each
        parameter.


    Examples
    --------
    >>> from sora.lightcurve.occdetect import occ_detect
    >>> params = occ_detect(flux, dflux, time, 0.2)
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

    # duration of the light curve
    time_span = time[-1] - time[0]
    if maximum_duration and (maximum_duration > time_span):
        warnings.warn('Occultation duration (maximum_duration={0}) '
                      'exceeds the time series length ({1:0.5f}).'
                      ' maximum_duration reset to the time series length.'
                      .format(maximum_duration, time_span))
        maximum_duration = time_span
    if not maximum_duration:
        maximum_duration = time_span

    if not dur_step:
        dur_step = cycle / 2

    if dur_step < cycle / 2:
        warnings.warn('The given dur_step is oversampled by a factor '
                      'of {0:0.1f} and has been reset to half a cycle '
                      '({1:0.4f}).'.format((cycle / 2.) / dur_step, cycle / 2.))
        dur_step = cycle / 2

    duration_grid = np.arange(dur_step, maximum_duration, dur_step)
    # initial occultation mask (all data points)
    mask = np.ones(len(time), dtype=bool)
    # initial detection rank
    rank = 1

    if snr_limit:
        # minimum SNR accepted in a detection for multiple search
        snr_value = snr_limit + 1
        occ0 = run_bls(flux, dflux, time, time_span, duration_grid)
        mask *= ~occ0['occ_mask']
        while snr_value > snr_limit:
            rank += 1
            occ1 = run_bls(flux, dflux, time, time_span, duration_grid, mask=mask,
                           rank=rank)
            if occ1['snr'] > snr_limit:
                snr_value = occ1['snr']
                mask *= ~occ1['occ_mask']
                occ0 = summarize_bls(occ0, occ1)
            else:
                snr_value = snr_limit

        if plot:
            plot_occ_detect(occ0, flux, time)
        return occ0
    elif n_detections:
        # search the n best fits
        occ0 = run_bls(flux, dflux, time, time_span, duration_grid)
        mask *= ~occ0['occ_mask']
        for i in range(n_detections - 1):
            rank += 1
            occ1 = run_bls(flux, dflux, time, time_span, duration_grid, mask=mask,
                           rank=rank)
            mask *= ~occ1['occ_mask']
            occ0 = summarize_bls(occ0, occ1)

        if plot:
            plot_occ_detect(occ0, flux, time)
        return occ0
    else:
        # search only the first best fit
        occ0 = run_bls(flux, dflux, time, time_span, duration_grid)

        if plot:
            plot_occ_detect(occ0, flux, time)

        return occ0


def plot_occ_detect(occ, flux, time):
    import matplotlib.pyplot as plt

    n = np.size(occ['rank'])
    if n > 1:
        # case for multiple detections
        plt.plot(time, flux, 'k.-')
        mask = np.zeros(len(time), dtype=bool)
        for i in range(n):
            trues = np.sum(occ['occ_mask'][i])
            plt.plot(time[occ['occ_mask'][i]], np.repeat(np.mean(flux[occ['occ_mask'][i]]), trues),
                     '.', label='Rank: ' + str(i + 1))
            mask += occ['occ_mask'][i]

        falses = list(mask).count(False)
        plt.plot(time[~mask], np.repeat(np.mean(flux[~mask]), falses), 'r.', label='Baseline')
        plt.xlabel('Time [seconds]')
        plt.ylabel('Relative Flux')
        plt.legend()

    else:
        # case for single occultation
        trues = list(occ['occ_mask']).count(True)
        falses = list(occ['occ_mask']).count(False)
        plt.plot(time, flux, 'k.-')
        plt.plot(time[occ['occ_mask']], np.repeat(np.mean(flux[occ['occ_mask']]), trues),
                 '.', label='Occultation')
        plt.plot(time[~occ['occ_mask']], np.repeat(np.mean(flux[~occ['occ_mask']]), falses),
                 'r.', label='Baseline')
        plt.xlabel('Time [seconds]')
        plt.ylabel('Relative Flux')
        plt.legend()


def run_bls(flux, dflux, time, per_grid, dur_grid, mask=None, rank=None):
    """ Private function to find the best box fit suitable to the data
    """
    from astropy.timeseries import BoxLeastSquares

    # object with no occultation mask
    mskmodel = BoxLeastSquares(time, flux, dy=dflux)
    # if there is no dflux array, reset it to None in case of
    # using a mask
    if dflux is None:
        dfluxmask = None
    else:
        dfluxmask = dflux[mask]

    # object with occultation mask
    if np.sum(mask):
        model = BoxLeastSquares(time[mask], flux[mask],
                                dy=dfluxmask)
    else:
        model = mskmodel

    r = model.power(per_grid, dur_grid, objective='snr', method='fast')
    # statistics of the BLS fit
    stats = model.compute_stats(r.period, r.duration, r.transit_time)
    # occultation mask of the event with respect to all data
    occ_mask = mskmodel.transit_mask(time, r.period, r.duration,
                                     r.transit_time)
    # parameters computation for clarity purposes
    occultation_duration = r.duration[0]
    central_time = stats['transit_times'][0]
    immersion_time = stats['transit_times'][0] - r.duration[0] / 2
    emersion_time = stats['transit_times'][0] + r.duration[0] / 2
    time_err = np.median(time[1:-1] - time[0:-2]) / 2
    depth = np.mean(flux[~occ_mask]) - np.mean(flux[occ_mask])
    depth_err = np.std(flux[occ_mask], ddof=1)
    baseline = np.mean(flux[~occ_mask])
    baseline_err = np.std(flux[~occ_mask], ddof=1)
    # If there is only one measurement during the occultation it will
    # use the baseline_err to compute SNR, otherwise it will use depth_err
    if np.sum(occ_mask) < 2:
        snr = depth / baseline_err
    else:
        snr = depth / depth_err

    # define rank
    if rank:
        rank = rank
    else:
        rank = 1

    return {'rank': rank,
            'occultation_duration': occultation_duration,
            'central_time': central_time,
            'immersion_time': immersion_time,
            'emersion_time': emersion_time,
            'time_err': time_err,
            'depth': depth, 'depth_err': depth_err,
            'baseline': baseline, 'baseline_err': baseline_err,
            'snr': snr, 'occ_mask': occ_mask}


def summarize_bls(dict1, dict2):
    """ Function to merge dictionaries returned by BLS and
        to keep values of common keys in list.
    """
    dict3 = {}
    for key, value in dict1.items():
        if key == 'occ_mask':
            sz = int(np.size(dict1[key]) / np.size(dict2[key]))
            if sz > 1:
                k = [None] * (sz + 1)
                for i in range(sz):
                    k[i] = dict1[key][i]
                k[i + 1] = dict2[key]
                dict3[key] = k
            else:
                dict3[key] = [dict1[key], dict2[key]]
        else:
            dict3[key] = np.append(dict1[key], dict2[key])
    return dict3
