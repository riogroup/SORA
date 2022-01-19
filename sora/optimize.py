from scipy import optimize
from scipy.integrate import quad
from numpy import inf, dot, linalg, sqrt, array, diag, ones, std, arange, isfinite, isinf, shape,sum
from numpy.random import choice
from warnings import warn
from copy import deepcopy
from tqdm import tqdm
from tqdm.notebook import tqdm_notebook


class OptimizeResult(dict):
    '''Descritpion
    
    Attributes
    ----------
    *params : Parameters object
        The best-fit parameters computed from the fit.
    method : str
        Method used to compute the best fit solution.
    var_names : list
        Ordered list of variable parameter names used in optimization, and
        useful for understanding the values in :attr:`initial_values` and
        :attr:`covar`.
    covar : numpy.ndarray
        Scaled covariance matrix from minimization, with rows and columns
        corresponding to :attr:`var_names`.
    initial_values : numpy.ndarray
        Dictionary of initial values for varying parameters.
    best_fit : numpy.ndarray
        List containing best fit values corresponding to :attr:`var_names`.
    std : numpy.ndarray
        Array containing estimated 1-sigma uncertainties, otherwise inf. 
        Corresponding to :attr:`var_names`.
    bootstrap: None
        Array containing the bootstraped solutions for each variable (col), by defalut None.
    success : str
        True if optimization algorithm converged.
    nvars : int
        Number of free variables used in fit.
    ndata : int
        Number of data points.
    dof : int
        Degrees of freedom in fit (ndata - nvars).
    residual : numpy.ndarray
        Residual array :math:`{\rm Resid_i}`. Return value of the objective
        function when using the best-fit values of the parameters.
    chisqr : float
         A weighted sum of squared deviations. If the objective function does
         not provide weighed residuals then it only expresses the sum of squared deviations.
    redchi : float
        The Reduced chi-square value: i is defined as chi-square per degree of freedom.
    lib : describe
    
    *emcee :  describe
    '''
    
    def __init__(self, **kwargs):
        for key, value in kwargs.items():
            setattr(self, key, value)

    def summary(self):
        s = ''
        print('FIT RESULTS:')
        for key in self.params.keys():
            ss = ''
            ss += '{: <20s}'.format(key[0:len(self.params[key].name) if (len(self.params[key].name) < 20) else 19])
            col = shape(self.params[key].std)

            if not self.params[key].free:
                    ss += 'value: {:<6.4g}  (Fixed)     '.format(self.params[key].value)
            else:
                if (len(col) > 0) and (col[0] > 1):
                    ss += 'value: {:<6g} (-{:.3g}, +{:.3g}) '.format(self.params[key].value, self.params[key].std[0], self.params[key].std[1])
                if (len(col) == 0):
                    ss += 'value: {:<6g} (+/-{:>.3g}) '.format(self.params[key].value, self.params[key].std)

                ss += ' (bounds [{:>g}, {:>g}]) '.format(self.params[key].minval, self.params[key].maxval)
            print(ss)

        print('\nFIT STATISTICS:')
        print(f'Method:           {self.method}')
        print(f'Data points:      {self.ndata}')
        print(f'Fitted variables: {self.nvars}')
        print(f'Deg. of freedom:  {self.dof}')
        print(f'Chi squared:      {self.chisqr:6g}')
        print(f'Red. Chi squared: {self.redchi:6g}')
        print(f'Success:          {self.success}')
        # print(f'Cov Matrix:       {"True" if hasattr(self,'covar') and self.covar is not None else "False"}')
        print(f'Bootstrap CI:     {"True" if self.bootstrap is not None else "False"}')
        return None
            

    def __repr__(self):
        ''''''
        # not returning string fix later
        for key, value in self.__dict__.items():
            print('{:<10s}: {}'.format(key, value))
        return ''
    
    def __str__(self):
        self.__repr__()
        return ''


    
def _fit_results(result, parameters, residual, bootstrap=None, lib='scipy', sigma=1, method=None):
    '''s.append(' value: {:.6g} '.format(self.params[key].value)
    '''

    # generate outputs
    output = OptimizeResult()

    # var_names
    output.var_names = parameters.get_names()

    # initial_values
    output.initial_values = array(parameters.get_varys())

    # sucess
    if hasattr(result, 'success'):
        output.success = result.success
    else:
        output.sucess = False

    # best_fit
    output.best_fit = result.x

    # nvars
    output.nvars = len(result.x)

    # ndata
    output.ndata = len(residual)

    # dof
    output.dof = len(residual) - len(result.x)

    # residual
    output.residual = array(residual)

    # chisq
    output.chisqr = array(residual).sum()

    # redchi
    output.redchi = array(residual).sum()/float(output.dof)

    # set uncertainties as inf by default
    # if other uncertainties are computed they will be replaced
    stderr = ones(len(result.x))*inf

    output.method = method
    # covar_matrix
    if hasattr(result, 'jac') and (method == 'least squares'):
        # try compute cov matrix
        try:
            # get approx Hessian by 2*Jac^T x Jac
            hessian = dot(result.jac.T, result.jac)
            # invert approx Hessian matrix
            inv_hessian = linalg.inv(hessian)
            output.covar = inv_hessian*output.redchi
            stderr = sqrt(diag(output.covar))*sigma
        except:
            output.covar = None
            warn(f'It was not possible to compute the covariance matrix.')
            pass

    # bootstrap
    output.bootstrap = bootstrap
   

    if bootstrap is not None:
        one_tailed, _ = quad(lambda x : (1.0 / np.sqrt(2*np.pi))*np.exp((-x**2) / 2.0), np.NINF, sigma)
        quantile = one_tailed-0.5
        low = np.quantile(bootstrap, 0.5-quantile, axis=0)
        high = np.quantile(bootstrap, 0.5+quantile, axis=0)
        stderr = []
        for i in range(len(low)):
            stderr.append([result.x[i]-low[i], high[i]-result.x[i]])
        stderr = np.array(stderr)
    
    # uncertainties
    output.std = stderr

    # scipy
    if lib == 'scipy':
        output.scipy = result

    # emcee
    # if method == 'emcee':
        # output.emcee = result

    # params
    var_index = 0
    params = deepcopy(parameters)
    for key in params.keys():
        params[key].initial_value = params[key].value
        params[key].std = inf
        if params[key].free:
            params[key].value = output.best_fit[var_index]
            params[key].std = stderr[var_index]
            var_index += 1

    output.params = params
    return output

    
    
    
def _bypass_fixed_variables(vary, func, parameters, *args, **kwargs):
        vary_index = 0
        params = deepcopy(parameters)
        for key, _ in params.items():
            if params[key].free:
                params[key].value = vary[vary_index]
                vary_index += 1

        return func(params, *args, **kwargs)


    
def _refactor_func_args(func, parameters, *args, index=None):
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
    
    accepted_methods = ['lm', 'trf', 'dogbox', 'differential_evolution']
    
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
                    raise ValueError(f'One or more values in the initial varying parameters are outside the provided bounds.')
    
    # case for LM
    if (method == 'lm'):
        bounds = (-inf, inf)
        
    # case for DIFFERENTIAL EVOLUTION
    if (method == 'differential_evolution'):
        if (bounds is None):
            bounds = parameters.get_bounds(transposed=False)
        else:
            for i, v in enumerate(vary):
                if not (bounds[i][0] <= v <= bounds[i][1]):
                    raise ValueError(f'One or more values in the initial varying parameters are outside the provided bounds.')
  
    return vary, bounds
            



def least_squares(func, parameters, bounds=None, args=(), method='lm', bootstrap=None, sigma=1, **kwargs):
    '''
    least_squares [summary]

    [extended_summary]

    Parameters
    ----------
    func : function
        Objective function
    parameters : Parameters object dict
        Object containing the parameters to be fitted.
    args : tuple, optional
        Additional arguments to be passed to the function, by default ()
    method : str, optional
        Optimization method used to perform the fit, by default 'lm'
        'lm'    : Levenberg-Marquardt (unbounded)
        'trf'   : Trust Region Reflective algorithm (bounds required)
        'dogbox : dogleg algorithm with rectangular trust regions, typical
                  use case is small problems with bounds (bounds required)
    **kwargs : dict, optional
        Aditional optios to pass to :scipydoc:`optimize.least_squares`

    '''

    
    # collect bootstrap variable status and pop it to avoid conflict with scipy
    if hasattr(kwargs, 'bootstrap'):
        bootstrap = kwargs['bootstrap']
        kwargs.pop('bootstrap')

    
    # get initial values and check external bonds provided
    vary, bounds = _prepare_fit(parameters, method, bounds)
    
   
    # check if initial parameters results are finite
    if np.isfinite(_bypass_fixed_variables(vary, func, parameters, *args)).sum() == 0:
        raise ValueError(f'Residuals are not finite at initial point. Check your initial parameters and/or objective function.')

    
    refact_args = _refactor_func_args(func, parameters, *args)

    # try execute the fit
    try:
        solution = optimize.least_squares(_bypass_fixed_variables, vary, args=refact_args, bounds=bounds, method=method, **kwargs)   
        residual = _bypass_fixed_variables(solution.x, func, parameters, *args)
        
        # if set compute bootstrap statistics
        bootstrap_array = None
        if isinstance(bootstrap, (float, int)):
            bootstrap_array = np.zeros((int(bootstrap), len(vary)))
            index = arange(len(args[0]))
            for boot_index in tqdm(range(int(bootstrap)), desc='Bootstraping'):
                new_index = choice(index, size=len(args[0]), replace=True, p=None)
                refact_args = _refactor_func_args(func, parameters, *args, index=new_index)
                # try to perform the fit
                try:
                    solution_bootstrap = optimize.least_squares(_bypass_fixed_variables, vary, args=refact_args, bounds=bounds, method=method, **kwargs)
                    bootstrap_array[boot_index,:] = solution_bootstrap.x
                except:
                    raise ValueError(f'An error occured during boostraping solutions.')
                
        return _fit_results(solution, parameters, residual, bootstrap=bootstrap_array, lib='scipy', method='least squares', sigma=sigma)
    except:
        raise ValueError(f'The optimization procedure failed.')
        


        
        
        

        
def differential_evolution(func, parameters, bounds=None, args=(), bootstrap=None, sigma=1, **kwargs):

        
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
        # if set compute bootstrap statistics
        bootstrap_array = None
        if isinstance(bootstrap, (float, int)):
            bootstrap_array = np.zeros((int(bootstrap), len(vary)))
            index = arange(len(args[0]))
            for boot_index in tqdm(range(int(bootstrap)), desc='Bootstraping'):
                new_index = choice(index, size=len(args[0]), replace=True, p=None)
                refact_args = _refactor_func_args(func, parameters, *args, index=new_index)
                # try to perform the fit
                try:
                    solution_bootstrap = optimize.differential_evolution(lambda *args: np.sqrt(np.sum(_bypass_fixed_variables(*args)**2)), bounds, args=refact_args, **kwargs)   
                    bootstrap_array[boot_index,:] = solution_bootstrap.x
                except:
                    raise ValueError(f'An error occured during boostraping solutions.')
                
        return _fit_results(solution, parameters, residual, bootstrap=bootstrap_array, method='differential_evolution', lib='scipy', sigma=sigma)
    except:
        raise ValueError(f'The optimization procedure failed.')       