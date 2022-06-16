from scipy import optimize
from scipy.integrate import quad
import numpy as np
from warnings import warn
from copy import deepcopy
from .core import Parameters


class NoTqdm():
    '''A dummy class to handle pbar when not showing progress bar.'''
    def __init__(self, *args, **kwargs):
        pass

    def update(self, *args, **kwargs):
        pass

    def close(self, *args, **kwargs):
        pass


class OptimizeResult(dict):
    '''An object that contains the results obtained by the fitting procedure
    
    Attributes
    ----------
    params : `Parameters` object
        The best-fit parameters computed from the fit.
    method : `str`
        Method used to compute the best fit solution.
    var_names : `list`
        Ordered list of variable parameter names used in optimization, and
        useful for understanding the values in :attr:`initial_values` and
        :attr:`covar`.
    covar : `numpy.ndarray`
        Scaled covariance matrix from minimization, with rows and columns
        corresponding to :attr:`var_names`. The covariance matrix is obtained 
        from the approximation of the inverse of the Hessian matrix. 
        This covariance matrix is scaled by the reduced chisqr.
    initial_values : `numpy.ndarray`
        Dictionary of initial values for varying parameters.
    best_fit : `numpy.ndarray`
        List containing best fit values corresponding to :attr:`var_names`.
    std : `numpy.ndarray`
        Array containing estimated N-sigma uncertainties, otherwise `inf`. 
        Corresponding to :attr:`var_names`. In the `least_squares` case the 
        formal numeric uncertainties are then derived from the covariance matrix.
        When uncertainties are unavailable Bootstraping can be used to derive the
        uncertainties. In this case confidence intervals (CI) limited by the provided 
        sigma (z-score), e.g., 1-sigma CI will return the [0.16,0.84] quantiles as 
        the low and high confidence intervals value.
    bootstrap: `None`, `numpy.ndarray`
        If not None, an array containing the bootstraped solutions for each variable (col), by defalut None.
    success : `str`
        True if optimization algorithm converged.
    nvars : `int`
        Number of free variables used in fit.
    ndata : `int`
        Number of data points.
    dof : `int`
        Degrees of freedom in fit (ndata - nvars).
    residual : `numpy.ndarray`
        Residual array of the objective function when using the parameters best-fit solution.
    chisqr : `float`
         A weighted sum of squared deviations. If the objective function does
         not provide weighed residuals then it only expresses the sum of squared deviations.
    redchi : float
        The Reduced chi-square value: i is defined as chi-square per degree of freedom.
    scipy : `scipy.OptimizeResult` object
        When scipy module is used the `scipy.OptimizeResult` object is also returned.
    emcee : `describe`
        Describe

    Returns
    -------
    object : `OptimizeResult`
        A OptimizeResult object containing the collection of the results obtained by the fit.
        '''
    
    def __init__(self, **kwargs):
        '''
        Contructor method of the OptimizeResults object.
        '''
        for key, value in kwargs.items():
            setattr(self, key, value)

    
    # THIS METHOD IS NOT RETURNING A STRING AND NEEDS A BETTER SOLUTION
    def __repr__(self):
        ''''''
        # not returning string fix later
        for key, value in self.__dict__.items():
            print('{:<10s}: {}'.format(key, value))
        return ''
    
    
    # THIS METHOD IS NOT RETURNING A STRING AND NEEDS A BETTER SOLUTION
    def __str__(self):
        ''''''
        self.__repr__()
        return ''


    def summary(self):
        '''
        Prints a summary of the results obtained by the fit contained in the
        OptimizeResult object.
        '''

        s = ''
        print('FIT RESULTS:')
        for key in self.params.keys():
            ss = ''
            ss += '{: <20s}'.format(key[0:len(self.params[key].name) if (len(self.params[key].name) < 20) else 19])
            col = np.shape(self.params[key].std)

            if not self.params[key].free:
                    ss += 'value: {:<6.4g}  (Fixed)     '.format(self.params[key].value)
            else:
                if (len(col) > 0) and (col[0] > 1):
                    ss += 'value: {:<6g} (-{:.3g}, +{:.3g}) '.format(self.params[key].value, self.params[key].std[0], self.params[key].std[1])
                if (len(col) == 0):
                    ss += 'value: {:<6g} (+/-{:>.3g}) '.format(self.params[key].value, self.params[key].std)

            if (self.params[key]):
                ss += ' (bounds [{:>g}, {:>g}]) '.format(self.params[key].minval, self.params[key].maxval)
            print(ss)

        print('\nFIT STATISTICS:')
        print(f'Method:           {self.method}')
        print(f'Data points:      {self.ndata}')
        print(f'Fitted variables: {self.nvars}')
        print(f'Deg. of freedom:  {self.dof}')
        print(f'Chi squared:      {self.chisqr:6g}')
        print(f'Red. Chi squared: {self.redchi:6g}')
        # print(f'Success:          {self.success}')
        # print(f'Cov Matrix:       {"True" if hasattr(self,'covar') and self.covar is not None else "False"}')
        print(f'Bootstrap CI:     {"True" if self.bootstrap is not None else "False"}')
        return None
            

    

# The following lines contain the fitting functions
    
def _scipy_results(result, parameters, residual, bootstrap=None, lib='scipy', sigma=1, method='Unavailabe'):
    '''
    Retuns an OptimizeResult object with the results obtained by the scipy module fitting.

    Parameters
    ----------
    result : `scipy.OptimizeResult` object
        Object with the results obtained by the fitting library
    parameters : `Parameters` object
        The Parameters object used to produce the fit.
    residual : `numpy.ndarray`
        Array containing the residuals of the objective function.
    bootstrap : `None`, `numpy.ndarray`, optional
        [description], by default None
    lib : `str`, mandatory
        Library used to produce the fit, by default 'scipy'
    sigma : `float`, `int`, optional
        Number of standard deviations corresponding the returned 
        uncertainties, by default 1
    method : `str`, optional
        Fitting method, by default `Unavailable`
        Options: `least_squares` or `differential_evolution`

    Returns
    -------
    object : `OptimizeResult`
        OptimizeResult object with the results obtained by the fitting.
    '''
    
    # instantiate the OptimizeResult object
    output = OptimizeResult()
    # add var_names attribute
    output.var_names = parameters.get_names()
    # add initial_values attribute
    output.initial_values = np.array(parameters.get_varys())
    # add success status attribute
    if hasattr(result, 'success'):
        output.success = result.success
    else:
        output.sucess = False
    # add best_fit parameters value attribute
    output.best_fit = result.x
    # add number of free variables attribute
    output.nvars = len(result.x)
    # add number of data points attribute
    output.ndata = len(residual)
    # add degrees of freedom attribute
    output.dof = len(residual) - len(result.x)
    #  add residual array attribute
    output.residual = np.array(residual)
    # add chisqr attribute
    output.chisqr = np.array(residual).sum()
    # add redchi attribute
    output.redchi = np.array(residual).sum()/float(output.dof)
    # set uncertainties as inf by default
    # if other uncertainties are computed they will be replaced
    stderr = np.ones(len(result.x))*np.inf
    # add method attribute
    output.method = method
    # if exists, compute the covariance matrix
    if hasattr(result, 'jac') and (method == 'least_squares'):
        # try compute cov matrix
        try:
            # get approx Hessian by 2*Jac^T x Jac
            hessian = np.dot(result.jac.T, result.jac)
            # invert approx Hessian matrix
            inv_hessian = np.linalg.inv(hessian)
            output.covar = inv_hessian*output.redchi
            stderr = np.sqrt(np.diag(output.covar))*sigma #sigma is normalized in this case
        except:
            output.covar = None
            warn(f'It was not possible to compute the covariance matrix.')
            pass

    # add bootstrap attribute
    output.bootstrap = bootstrap
    # compute uncertainties from bootstrap
    if bootstrap is not None:
        one_tailed, _ = quad(lambda x : (1.0 / np.sqrt(2*np.pi))*np.exp((-x**2) / 2.0), np.NINF, sigma)
        quantile = one_tailed-0.5
        low = np.quantile(bootstrap, 0.5-quantile, axis=0)
        high = np.quantile(bootstrap, 0.5+quantile, axis=0)
        stderr = []
        for i in range(len(low)):
            stderr.append([result.x[i]-low[i], high[i]-result.x[i]])
        stderr = np.array(stderr)
    
    # add std attribute
    output.std = stderr
    # add scipy attribute
    if lib == 'scipy':
        output.scipy = result

    #FOR FUTURE IMPLEMENTATION
    # emcee
    # if method == 'emcee':
        # output.emcee = result

    # add params attribute and update results
    var_index = 0
    params = deepcopy(parameters)
    for key in params.keys():
        params[key].initial_value = params[key].value
        params[key].std = np.inf
        if params[key].free:
            params[key].value = output.best_fit[var_index]
            params[key].std = stderr[var_index]
            var_index += 1

    output.params = params
    
    return output


def _fastchi_results(result, parameters, residual, n_samples, sigma_range, sigma, sigma_samples, ndata):
    '''
    Retuns an OptimizeResult object with the results obtained by the scipy module fitting.

    Parameters
    ----------
    result : `numpy.ndarray`
        _marching_grid function result.
    parameters : `Parameters` object
        The Parameters object used to produce the fit.
    

    Returns
    -------
    object : `OptimizeResult`
        OptimizeResult object with the results obtained by the fitting.
    '''
    
    # instantiate the OptimizeResult object
    output = OptimizeResult()
    # add var_names attribute
    output.var_names = parameters.get_names()
    # add initial_values attribute
    output.initial_values = np.array(parameters.get_varys())
    
    samples_full, chisqr_full = result
    tsamples_full = samples_full.T

    best_fit_index = np.argmin(chisqr_full)
    samples_index = np.argwhere(chisqr_full <= (chisqr_full.min() + sigma**2)).T[0]
    samples_output_index = np.argwhere(chisqr_full <= (chisqr_full.min() + sigma_range**2)).T[0]
    
    # discard exceding samples when s
    if sigma_samples is None:            
        tmp_i = 0
        while (len(samples_index) > n_samples):
            if (tmp_i == best_fit_index):
                tmp_i += 1
            else:
                samples_index = np.delete(samples_index, tmp_i)
    
    samples, chisqr = samples_full[samples_index].T, chisqr_full[samples_index]
    samples_output, chisqr_output = samples_full[samples_output_index].T, chisqr_full[samples_output_index]
    best_fit_index = np.argmin(chisqr)
    
    # add best_fit parameters value attribute
    output.best_fit = samples[:,best_fit_index]
    # add number of free variables attribute
    output.nvars = len(samples[:,best_fit_index])
    # add number of data points attribute
    output.ndata = ndata
    # add degrees of freedom attribute
    output.dof = output.ndata - output.nvars
    #  add residual array attribute
    output.residual = np.array(residual)
    # add chisqr attribute
    output.chisqr = chisqr[best_fit_index]
    # add redchi attribute
    output.redchi = output.chisqr/output.dof
    # add bootstrap attribute
    output.bootstrap = None
    # add sample values attribute
    output.samples = samples_output
    # add chi values attribute
    output.chi = chisqr_output
    # if other uncertainties are computed they will be replaced
    stderr = []
    for n in range(output.nvars):
        tmp_idx = (chisqr < output.chisqr + sigma**2)
        left = min(samples[n][tmp_idx])
        tmp_idx = (chisqr < output.chisqr + sigma**2)
        right = max(samples[n][tmp_idx])
        stderr.append([abs(output.best_fit[n]-left), abs(right-output.best_fit[n])])

    # add method attribute
    output.method = 'Chisqr'
    # add sigma std level
    output.sigma = sigma
    # add sigma_range interval
    output.sigma_range = sigma_range
    # add std attribute
    output.std = stderr
    # add params attribute and update results
    var_index = 0
    params = deepcopy(parameters)
    for key in params.keys():
        params[key].initial_value = params[key].value
        params[key].std = np.inf
        if params[key].free:
            params[key].value = output.best_fit[var_index]
            params[key].std = stderr[var_index]
            var_index += 1

    output.params = params
    
    return output    
   
    
def _bypass_fixed_variables(vary, func, parameters, *args, **kwargs):
    '''
    Bypass non free variables during the fitting procedure.

    Parameters
    ----------
    vary : `np.ndarray`
        Array of free variables.
    func : `userfunc`
        User provided function.
    parameters : `Parameters` object
        Parameters used during the fitting

    Returns
    -------
    Bypassed `userfunc` with fixed values *args and **kwargs.
    '''
    vary_index = 0
    params = deepcopy(parameters)
    for key, _ in params.items():
        if params[key].free:
            params[key].value = vary[vary_index]
            vary_index += 1

    return func(params, *args, **kwargs)


    
def _refactor_func_args(func, parameters, *args, index=None):
    '''
    Refactor function arguments

    [extended_summary]

    Parameters
    ----------
    func : `userfunc`
        User provided function.
    parameters : `Parameters` object
        Parameters used during the fitting
    index : `None`, `list`, optional
        When not None provides the indices for bootstraping, by default None
    '''
        # refactor arguments
    refactored_args = []
    refactored_args.append(func)
    refactored_args.append(parameters)
    if index is not None:
        for value in args:
            refactored_args.append(value[index])
    else:
        for value in args:
            refactored_args.append(value)
    return refactored_args



def _prepare_fit(parameters, method, bounds):
    '''
    Prepares the free variables values and the checks for 
    the appropriate bounds used in each method.

    Parameters
    ----------
    parameters : `Parameters` object
        Parameters used during the fitting
    method : `str`
        The method used to produce the fit.
        Options: `lm`, `trf`, `dogbox`, and `differential evolution`.
    bounds : `tuple`, `scipy.optimize.Bounds` object
        If not None are the bounds externaly provided by the user.

    Returns
    -------
    [values, bounds] : `list`
        Values of the free variables and bounds in the appropriate format for
        each fitting method.

    Raises
    ------
    ValueError
        When `method`, `parameters` or `bounds` are incorrectly provided.
    '''
    
    accepted_methods = ['lm', 'trf', 'dogbox', 'differential_evolution', 'fastchi']
    
    # chec if the method is available in the routines
    if not any((method == a) for a in accepted_methods):
        raise ValueError(f'An invalid method `{method}` was provided. Please refer to the documentation.')
    
    # check if the parameters provided are an Parameters instantiated object
    if parameters is not None and not isinstance(parameters, Parameters):
        raise ValueError(f'{parameters} is not a valid Parameters object dictionary.')

    # get varying parameters
    vary = parameters.get_varys()

    # for each type of method ensure the correct bounds format 
    # case for TRF and DOGBOX
    if (method == 'trf') or (method == 'dogbox'):
        if (bounds is None):
            bounds = parameters.get_bounds(transposed=True)
        else:
            for i, v in enumerate(vary):
                if not (bounds[0][i] <= v <= bounds[1][i]):
                    raise ValueError(f'One or more values in the initial varying parameters is outside the provided bounds.')
    
    # case for LM
    if (method == 'lm'):
        bounds = (-np.inf, np.inf)

    # case for DIFFERENTIAL EVOLUTION
    if (method == 'differential_evolution') or (method == 'fastchi'):
        if (bounds is None):
            bounds = parameters.get_bounds(transposed=False)
        else:
            for i, v in enumerate(vary):
                if not (bounds[i][0] <= v <= bounds[i][1]):
                    raise ValueError(f'One or more values in the initial varying parameters are outside the provided bounds.')

    return vary, bounds
            


def least_squares(func, parameters, bounds=None, args=(), algorithm='lm', bootstrap=None, sigma=1, show_progress=False, **kwargs):
    '''
    Performs the optimization of function parameters using convervenge algorithms and the scipy module.

    Parameters
    ----------
    func : `userfun`
        Objective function provided by the user. The objective function returns the 
        squared (or weighed squared) residual array, e.g., (model - data )^2, and must
        take the arguments as follows: `userfunc(parameters, *args, *kwargs)`
    parameters : `Parameters` object dict
        Object containing the parameters to be fitted.
    bounds : `None`, `tuple`
        Option to provide an array containing bounds external to the Parameters
        object. For more information on how to create bounds arrays 
        see :scipydoc:`optimize.Bounds` and `optimize.least_squares`
    args : `tuple`, optional
        Arguments passed to the `userfun`, by default ()
    algorithm : `str`, optional
        Optimization algorithm used to perform the fit, by default 'lm'
        'lm'    : Levenberg-Marquardt (unbounded).
        'trf'   : Trust Region Reflective algorithm (bounds required).
        'dogbox : dogleg algorithm with rectangular trust regions, typical
                  use case is small problems with bounds (bounds required).
    bootstrap: `None`, `float`, `int`
        The number of bootstrap computations desired to compute the confidence
        intervals of the results, by default None.
    sigma : float, int
        The number of sigmas correponding to the error bars obtained from the
        covariance matrix or with bootstraping.
    show_progress : `bool`
        If `True` shows a progress bar relative to the boostraping calculations.
        Depends on `tqdm` module to function. By default, False. 
    **kwargs : dict, optional
        Aditional options to pass see :scipydoc:`optimize.least_squares`

    
    Returns
    -------
    object : `OptimizeResult`
        OptimizeResult object with the results obtained by the fitting.
    '''
    
    # collect bootstrap variable status and pop it to avoid conflict with scipy
    if hasattr(kwargs, 'bootstrap'):
        bootstrap = kwargs['bootstrap']
        kwargs.pop('bootstrap')

    # get initial values and check external bonds provided
    vary, bounds = _prepare_fit(parameters, algorithm, bounds)
   
    # check if initial parameters results are finite
    if np.isfinite(_bypass_fixed_variables(vary, func, parameters, *args)).sum() == 0:
        raise ValueError(f'Residuals are not finite at initial point. Check your initial parameters and/or objective function.')
    
    refact_args = _refactor_func_args(func, parameters, *args)

    # try execute the fit
    try:
        solution = optimize.least_squares(_bypass_fixed_variables, vary, args=refact_args, bounds=bounds, method=algorithm, **kwargs)   
        residual = _bypass_fixed_variables(solution.x, func, parameters, *args)
        
        # if set compute bootstrap statistics
        bootstrap_array = None
        if isinstance(bootstrap, (float, int)):
            bootstrap_array = np.zeros((int(bootstrap), len(vary)))
            index = np.arange(len(args[0]))
            
            # progress bar settings
            if show_progress:
                pbar = tqdm(total=bootstrap, desc='Bootstraping')
            else:
                pbar = NoTqdm()
                if show_progress:
                    warn(f'`tqdm` module is required to show progress bars.')

            # single thread processing
            for boot_index in range(int(bootstrap)):
                new_index = np.random.choice(index, size=len(args[0]), replace=True, p=None)
                refact_args = _refactor_func_args(func, parameters, *args, index=new_index)
                # try to perform the fit
                try:
                    solution_bootstrap = optimize.least_squares(_bypass_fixed_variables, vary, args=refact_args, bounds=bounds, method=algorithm, **kwargs)
                    bootstrap_array[boot_index,:] = solution_bootstrap.x
                except:
                    raise ValueError(f'An error occured during boostraping solutions.')
        
                pbar.update()

            pbar.close()
                
        return _scipy_results(solution, parameters, residual, bootstrap=bootstrap_array, lib='scipy', method='least_squares', sigma=sigma)
    except:
        raise ValueError(f'The optimization procedure failed.')
        
        
        
        

        
def differential_evolution(func, parameters, bounds=None, args=(), bootstrap=None, sigma=1, show_progress=False, **kwargs):
    '''
    Performs the optimization of function parameters using `differential_evolution` method from
    scipy optimization module.

    Parameters
    ----------
    func : `userfun`
        Objective function provided by the user. The objective function returns the 
        squared (or weighed squared) residual array, e.g., (model - data )^2, and must
        take the arguments as follows: `userfunc(parameters, *args, *kwargs)`
    parameters : `Parameters` object dict
        Object containing the parameters to be fitted.
    bounds : `None`, `tuple`
        Option to provide an array containing bounds external to the Parameters
        object. For more information on how to create bounds arrays 
        see :scipydoc:`optimize.Bounds` and `optimize.differential_evolution`
    args : `tuple`, optional
        Arguments passed to the `userfun`, by default ()
    bootstrap: `None`, `float`, `int`
        The number of bootstrap computations desired to compute the confidence
        intervals of the results, by default None.
    sigma : float, int
        The number of sigmas correponding to the error bars obtained from the
        covariance matrix or with bootstraping.
    show_progress : `bool`
        If `True` shows a progress bar relative to the boostraping calculations.
        Depends on `tqdm` module to function. By default, False. 
    **kwargs : dict, optional
        Aditional options to pass during the fitting see :scipydoc:`optimize.differential_evolution`

    
    Returns
    -------
    object : `OptimizeResult`
        OptimizeResult object with the results obtained by the fitting.
    '''
        
    # collect bootstrap variable status and pop it to avoid conflict with scipy
    if hasattr(kwargs, 'bootstrap'):
        bootstrap = kwargs['bootstrap']
        kwargs.pop('bootstrap')
    
    # get initial values and check external bonds provided
    vary, bounds = _prepare_fit(parameters, 'differential_evolution', bounds) 
   
    # check if initial parameters results are finite
    if np.isfinite(_bypass_fixed_variables(vary, func, parameters, *args)).sum() == 0:
        raise ValueError(f'Residuals are not finite at initial point. Check your initial parameters and/or objective function.')
    
    refact_args = _refactor_func_args(func, parameters, *args)
    # try execute the fit
    try:
        solution = optimize.differential_evolution( lambda *args: np.sqrt(np.sum(_bypass_fixed_variables(*args)**2)), bounds, args=refact_args, **kwargs)   
        residual = _bypass_fixed_variables(solution.x, func, parameters, *args)
        
        # if set compute bootstrap statistics
        bootstrap_array = None
        if isinstance(bootstrap, (float, int)):
            bootstrap_array = np.zeros((int(bootstrap), len(vary)))
            index = np.arange(len(args[0]))
            
            # progress bar settings
            if show_progress:
                pbar = tqdm(total=bootstrap, desc='Bootstraping')
            else:
                pbar = NoTqdm()
                if show_progress:
                    warn(f'`tqdm` module is required to show progress bars.')

            # single thread processing
            for boot_index in range(int(bootstrap)):
                new_index = np.random.choice(index, size=len(args[0]), replace=True, p=None)
                refact_args = _refactor_func_args(func, parameters, *args, index=new_index)
                # try to perform the fit
                try:
                    solution_bootstrap = optimize.differential_evolution(lambda *args: np.sqrt(np.sum(_bypass_fixed_variables(*args)**2)), bounds, args=refact_args, **kwargs)   
                    bootstrap_array[boot_index,:] = solution_bootstrap.x
                except:
                    raise ValueError(f'An error occured during boostraping solutions.')
                
                pbar.update()
            
            pbar.close()
                
        return _scipy_results(solution, parameters, residual, bootstrap=bootstrap_array, method='differential_evolution', lib='scipy', sigma=sigma)
    except:
        raise ValueError(f'The optimization procedure failed.')       


def _marching_grid(func, initpars, bounds, args=(), sigma_range=3, sigma=1, sigma_samples=None, samples=500, marching_grid=True, run_size=10000, threads=None, show_progress=False):
    '''
    Multidimensional marching grid algorithm to optimize brute force minimum regions finding.

    Parameters
    ----------
    func : `userfunc`
        Objective function used to return the resiudals.
    bounds : `tuple`, `list`
        An array containing the bounding limits for each parameter. It should
        have the shape (Npar,2). Example: `bounds = ((-10,10), (0,1), (3,5))`.
    args : `tuple`, optional
        Aditional arguments passed to the user function, by default ().
    sigma_range: `int`, optional
        Interval range for sampling. All returned samples will lie
        within 0 and sigma_range standard deviations, by default 3.
    sigma : `int`, optional
        Number of standard deviations of the returned uncertaties, by default 1.
    sigma_samples : `int`, optional
        The minimum number of samples used to estimate the returned uncertainties
        defined by `sigma`, by default None.
        When defined it superseeds the `samples` variable since it will return
        most likely a higher number of samples.
    samples : `int`, optional
        Number of samples within the provided `sigma_range`, by default 2000.
    marching_grid : `bool`, optional
        If `True` it will use the marching grid approach to march faster towards the best solution.
        When `False` sampling is done uniformly within the boundings provided, by default True.
    run_size : `int`, optional
        Number of simulations performed within each sampling run. By default 10000.
    threads : `int`
        If multithreading is desired it sets the number of parallel processes, by default None.
    show_progress : `bool`
        If `True` shows a progress bar relative to the calculations. 
        Depends on `tqdm` module to function. By default, False. 

    Returns
    -------
    [parameters, residual] : `np.ndarray`
        A `numpy.ndarray` containing the sampled parameters and residuals.
    '''
    
    params, residual = [], []
    
    res = func(initpars,*args)
    params.append(list(initpars))
    residual.append(res.sum())

    # define the total iteration counter 
    # define the counter display and can be settet to be used with dependency check
    if show_progress:
        pbar = tqdm(total=samples, position=0)
    else:
        pbar = NoTqdm()
        if show_progress:
            warn(f'`tqdm` module is required to show progress bars.')
        
    if (sigma > sigma_range):
        sigma_range = sigma
        warn(f'`sigma_range` provided is below `sigma` requested, therefore set to the value of `sigma`.')
    
    counter = 0
    counter_previous = 0    
    sigma_counter = 0
    new_bounds = list(bounds)
    sigma_samples = 1 if sigma_samples is None else sigma_samples
    while ((counter < samples) or (sigma_counter < sigma_samples)):     
        # generate random samples
        p = []
        for i, b in enumerate(new_bounds):
            if (b[0] < bounds[i][0]):
                b[0] = bounds[i][0] 
            if (b[1] > bounds[i][1]):
                b[1] = bounds[i][1]
            
            p.append(np.random.uniform(low=b[0], high=b[1], size=run_size))
        p = np.array(list(p)).T

    
        # compute residuals for candidates multithreading
        if (threads is not None) and isinstance(threads, (float, int)):
            # create arguments pack
            pool_args = [ (pars, *args) for pars in p ]


            with Pool(processes=int(threads)) as pool:
                pool_res = pool.starmap(func, pool_args)

            for i in range(int(run_size)):
                params.append(list(p[i]))
                residual.append(pool_res[i].sum())
        
        # or instead compute residuals for candidates linearly
        else:
            if not isinstance(threads, (float, int)) and (threads is not None):
                warn(f'Provide the number of workers/threads to run in parallel in the keyword `threads`.')
            
            for i in range(int(run_size)):
                res = func(p[i],*args)
                params.append(list(p[i]))
                residual.append(res.sum())

        sigma_threshold = min(residual) + sigma_range**2
        min_res_index = ( residual <= sigma_threshold ) 

        #if there are at least two data points redefine bounds
        if marching_grid:
            idx_min = ( residual == min(residual) )
            params_check = np.array(params)
            for i, boundval in enumerate(bounds):
                leftlim =  min(params_check[min_res_index,i])
                rightlim = max(params_check[min_res_index,i])
                value = params_check[idx_min,i]              
                dv = max(value[0]-leftlim, rightlim-value[0])
                if (dv == 0):
                    new_bounds[i] = [boundval[0], boundval[1]]
                else:
                    new_bounds[i] = [leftlim-dv, rightlim+dv]


        sigma_counter = np.sum( residual <= (min(residual) + sigma**2))
        counter = np.sum( residual <= sigma_threshold )

        if (sigma_samples > 1) or (sigma_counter < sigma_samples):
            counter = np.round(samples*(sigma_counter/sigma_samples))

        pbar.update(counter-counter_previous if counter < samples else samples-counter_previous)

        counter_previous = counter

    pbar.close()
    return np.array(params), np.array(residual)


def fastchi(func, parameters, bounds=None, args=(), **kwargs):
    '''
    Performs the optimization of function parameters using `differential_evolution` method from
    scipy optimization module.

    Parameters
    ----------
    func : `userfun`
        Objective function provided by the user. The objective function returns the 
        squared (or weighed squared) residual array, e.g., (model - data )^2, and must
        take the arguments as follows: `userfunc(parameters, *args, *kwargs)`
    parameters : `Parameters` object dict
        Object containing the parameters to be fitted.
    bounds : `None`, `tuple`
        Option to provide an array containing bounds external to the Parameters
        object. For more information on how to create bounds arrays 
        see :scipydoc:`optimize.Bounds` and `optimize.differential_evolution`
    args : `tuple`, optional
        Arguments passed to the `userfun`, by default ()

    **kwargs :         
        Aditional options to pass during the fitting.
    sigma : `int`, optional
        Confidence limits within with the samples should lie, by default 1.    
    sigma_range: `int`, optional
        Confidence interval range for sampling. All returned samples will lie
        within 0 and sigma_range standard deviations, by default 3.
    samples : `int`, optional
        Number of simulations to be returned within the provided sigma, by default 2000.
    sigma_samples : `int`
        The minimum number of samples that lie within the defined `sigma` interval,
        by default 1.
    marching_grid : `bool`, optional
        If `True` it will use the marching grid approach to march faster towards the best solution.
        When `False` sampling is done uniformly within the boundings provided, by default True.
    run_size : `int`, optional
        Number of simulations performed within each sampling run. By default 10000.
    threads : `int`
        If multithreading is desired it sets the number of parallel processes, by default None.
    show_progress : `bool`
        If `True` shows a progress bar relative to the calculations. 
        Depends on `tqdm` module to function. By default, False. 

    Returns
    -------
    object : `OptimizeResult`
        OptimizeResult object with the results obtained by the fitting.
    '''
        
    
    # get initial values and check external bonds provided
    vary, bounds = _prepare_fit(parameters, 'fastchi', bounds) 
    

    # check if initial parameters results are finite
    if np.isfinite(_bypass_fixed_variables(vary, func, parameters, *args)).sum() == 0:
        raise ValueError(f'Residuals are not finite at initial point. Check your initial parameters and/or objective function.')
    
    refact_args = _refactor_func_args(func, parameters, *args)
    # try execute the fit
    try:
        solution = _marching_grid(_bypass_fixed_variables, vary, bounds, args=refact_args, **kwargs)
        tmp_samples, tmp_chisqr = solution   
        residual = _bypass_fixed_variables(tmp_samples[np.argmin(tmp_chisqr)], func, parameters, *args)
        #return solution
        return _fastchi_results(solution, parameters, residual, 
                                kwargs['samples'] if 'samples' in kwargs else 2000, 
                                kwargs['sigma_range'] if 'sigma_range' in kwargs else 3,
                                kwargs['sigma'] if 'sigma' in kwargs else 1, 
                                kwargs['sigma_samples'] if 'sigma_samples' in kwargs else None, len(residual))

    except:
        raise ValueError(f'The optimization procedure failed.')       