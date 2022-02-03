import numpy as np

# try import ipython display for table display in notebooks
try:
    from IPython.display import HTML, display
    from IPython import __version__
    HAS_IPYTHON = int(__version__[0]) >= 7
except ImportError:
    HAS_IPYTHON = False


def _in_ipynb():
    '''Check if the code is running on an ipython notebook'''
    try:
        cfg = get_ipython().config 
        if cfg['IPKernelApp']['parent_appname'] == 'ipython-notebook':
            return True
        else:
            return False
    except NameError:
        return False


class Parameter:
    '''A collection of values in an object to be used in a fitting procedure.
       The structure of this class heavely borrows its structure from parameters.py (lmfit)
    '''

    def __init__(self, name, value=None, minval=-np.inf, maxval=np.inf, free=True, std=None, initial_value=None):
        '''
        Constructor method for the Parameter object.

        Parameters
        ----------
        name : `str`
            A name or label to describe the parameter
        value : `float` 
            When first created it contains the initial estimate of the parameter. 
            In a result object it contains the best fit value, by default 0.
        minval : `float`, optional
            The lower bound used in the parameter variation, by default -inf
        maxval : `float`, optional
            The upper bound used in the parameter variation, by default inf
        free : bool, optional
            Defines if the parameter is allowed to vary, by default True
        std : `float`, `array`
            In the result object it contains the computed uncertainty.
        initial_value : `float`
            In the result object it contains the estimate before the fit.

        
        Returns
        -------
        object : `Parameter`
            A Parameter object containing the collection of values of a single parameter.

        '''

        self.name = name
        self.free = free
        self.std = std
        self.initial_value = initial_value
        self.minval = minval
        self.maxval = maxval

        if self.minval > self.maxval:
            self.minval, self.maxval = maxval, minval

        if np.isclose( self.minval, self.maxval, atol=1e-15, rtol=1e-15):
            raise ValueError(f'Parameter {self.name} has `minval` equal to maxval')
        
        if self.minval <= value <= self.maxval:
            self.value = value
        else:
            raise ValueError(f'Value `{value}` at parameter `{self.name}` is outside the provided bounds `({self.minval}, {self.maxval})`')
            
            
    def __repr__(self):
        """Return printable representation of a Parameter object."""
        s = []
        if not self.free:
            sval = f"value={self.value:.5g} (fixed)"
            s.append(sval)
        else:
            sval = f"value={self.value:.5g}"
            if self.std is not None:
                if np.size(self.std) == 2:
                    sval += f" -{self.std[0]:.3g}/+{self.std[1]:.3g}"
                if np.size(self.std) == 1:    
                    if np.isinf(self.std):
                        sval += f"-"
                    else:
                        sval += f" +/- {self.std:.3g}"
            s.append(sval)
            s.append(f", bounds=[{repr(self.minval)}:{repr(self.maxval)}]")
            s.append(f", free={repr(self.free)}")
        return f"{''.join(s)}"


    def update(self, **kwargs):
        for key, value in kwargs.items():
            setattr(self, key, value)





class Parameters(dict):
    '''
    Creates a dictionary of Parameter objects.
    The structure of this class heavely borrows its structure from parameters.py (lmfit)

    The dictionary contains all the parameters necessary to perform the fit of
    a model/objective function.

    Parameters
    ----------
    dict : `dict`
        A dictionary of Parameter objects.
    '''

    def __init__(self):
        '''
        Contructor method
        '''
        super().__init__(self)



    def _addpar(self, name, parameter):
        '''
        Private method to add the Parameter object to the Parameters dictionary.

        Parameters
        ----------
        name : `string`
            A name or label to describe the parameter
        parameter : `Parameter` object
            Parameter object containing the parameter values

        Raises
        ------
        ValueError
            Error if passed object is not a valid Parameter object.
        '''

        # check if the provided object is a valid instance of the Parameter object
        if parameter is not None and not isinstance(parameter, Parameter):
            raise ValueError(f'{parameter} is not a valid Parameter object.')
        # set item into dictionary
        dict.__setitem__(self, name, parameter)


    def add(self, name, value=None, minval=-np.inf, maxval=np.inf, free=True, std=None, initial_value=None):
        '''
        Include a Parameter values collection into the Parameters dictionary.

        Parameters
        ----------
        name : `str`
            A name or label to describe the parameter
        value : `float`, optional
            When first created it contains the initial estimate of the parameter. 
            In a result object it contains the best fit value, by default None.
        minval : `float`, optional
            The lower bound used in the parameter variation, by default -inf
        maxval : `float`, optional
            The upper bound used in the parameter variation, by default inf
        free : bool, optional
            Defines if the parameter is allowed to vary, by default True
        std : `float`, `array`
            In the result object it contains the computed uncertainty.
        initial_value : `float`
            In the result object it contains the estimate before the fit.

        
        Returns
        -------
        object : `Parameters`
            A Parameters object containing the collection of parameters.

        '''

        self._addpar(name, Parameter(name, value=value, minval=minval, maxval=maxval, free=free, std=std, initial_value=value) )


    def add_many(self, *parlist):
        '''
        Add many parameters using a list of tuples.

        Parameters
        ----------
        *parlist : :obj:`list` of :obj:`tuple`
            A list of tuples containing the parameters in a sequence.
            The order in each tuple must follow the sequence:
            ``(name, value, minval, maxval, free)``
        
        
        Returns
        -------
        object : `Parameters`
            A Parameters object containing the collection of parameters.

        '''

        for par in parlist:
            if not isinstance(par, Parameter):
                parobj = Parameter(*par)
                parobj.initial_value = parobj.value            
            self._addpar(parobj.name, parobj )



    def remove(self, name):
        '''
        Method to remove a Parameter from the Parameters dictionary.

        Parameters
        ----------
        name : `str`
            The name of the parameter to be removed.
        '''
        dict.pop(self, name)


    def get_bounds(self, transposed=False):
        '''
        Method to get the bounds recorded in the Parameters object.

        Parameters
        ----------
        transposed : bool, optional
            Returns the bounds in a transposed array, by default False.

        Returns
        -------
        array : `tuple`
            Tuple array containing the bounds values.
        '''
        
        if not transposed:
            bounds = []
            for key in self.valuesdict().keys():
                if self[key].free:
                    bounds.append((self[key].minval, self[key].maxval))
            return tuple(bounds)
        
        else:
            lowerbound = []
            upperbound = []
            for key in self.valuesdict().keys():
                if self[key].free:
                    lowerbound.append(self[key].minval)
                    upperbound.append(self[key].maxval)
            return (tuple(lowerbound),tuple(upperbound))


    def get_values(self):
        '''
        Method to get the values stored in the Parameters object.

        Returns
        -------
        list : `float`
            List of parameters values.
        '''

        values = []
        for key in self.valuesdict().keys():
            values.append(self[key].value)
        return values
    
    
    def get_varys(self):
        '''
        Method to get the FREE values stored in the Parameters object.

        Returns
        -------
        list : `tuple`
            List of FREE parameters values.
        '''
        values = []
        for key in self.valuesdict().keys():
            if self[key].free:
                values.append(self[key].value)
        return tuple(values)
    
    
    def get_uncertainties(self):       
        '''
        Method to get the uncertainties stored in the Parameters object.

        Returns
        -------
        list : `float`, `list`
            List of parameters uncertainties.
        '''
        values = []
        for key in self.valuesdict().keys():
            values.append(self[key].std)
        return values

    
    def get_names(self):
        '''
        Method to get the names of each parameter in the object.

        Returns
        -------
        list : `str`, `list`
            List of parameters names.
        '''
        values = []
        for key in self.valuesdict().keys():
            if self[key].free:
                values.append(self[key].name)
        return values
    
    
    def valuesdict(self):
        '''
        Method returns a dictionary of parameter and values.

        Returns
        -------
        dict : `dict`
            Dictionary of parameters names and values.
        '''
        return {p.name: p.value for p in self.values()}


    def __repr__(self):
        """__repr__ from OrderedDict."""
        
        # __repr__ if in a jupyter notebook
        if HAS_IPYTHON and _in_ipynb():

            for name in self.valuesdict():
                if self[name].std is None:
                    header = ['Parameter','Value','Bounds','Initial Value', 'Free']
                else:
                    header = ['Parameter','Value','Standard Error', 'Bounds','Initial Value', 'Free']

            s = []
            s.append('<table>')
            s.append('<tbody>')
            # header
            s.append('<tr>')
            for i, h in enumerate(header):
                if i == 0:
                    s.append('<td style="text-align: left;"><strong>'+str(h)+'</strong></td>')
                else:
                    s.append('<td><strong>'+str(h)+'</strong></td>')
            s.append('</tr>')
            # values
            for name in self.valuesdict():
                
                if not self[name].free:
                    s.append('<tr>')
                    s.append('<td style="text-align: left;">{} (fixed)</td>'.format(name))
                    s.append('<td>{:.5g}</td>'.format(self[name].value))
                    if self[name].std is not None:
                        s.append('<td>-</td><td>-</td><td>-</td><td>False</td>')
                    else:
                        s.append('<td>-</td><td>-</td><td>False</td>')
                else:
                    s.append('<tr>')
                    s.append('<td style="text-align: left;">{}</td>'.format(name))
                    s.append('<td>{:.5g}</td>'.format(self[name].value))
                    if self[name].std is not None:
                        if np.size(self[name].std) == 2:
                            s.append('<td>-{:.3g}/+{:.3g}</td>'.format(self[name].std[0], self[name].std[1]))
                        if np.size(self[name].std) == 1:
                            if np.isinf(self[name].std):
                                s.append('<td>-</td>')
                            else:
                                s.append('<td>+/-{:.3g}</td>'.format(self[name].std))
                    s.append('<td>{:.3g} : {:.3g}</td>'.format(self[name].minval, self[name].maxval))
                    s.append('<td>{}</td>'.format(self[name].initial_value))
                    s.append('<td>{}</td>'.format(self[name].free))
                    s.append('</tr>')
            s.append('</tbody>')
            s.append('</table>')

            s = ''.join(str(_) for _ in s)
            display(HTML(s))
            return ''

        # __repr__ if not in a jupyter notebook
        if not _in_ipynb():
            s = []
            for item in self.items():
                s.append(item)
                s.append('\n')
            s = ''.join(str(_) for _ in s)
            return f'{s}'        


    def __str__(self):
        """Return printable representation of a Parameters."""
        s = ['[ Parameters:\n']
        for item in self.items():
            s.append(item)
            s.append('\n')
        s.append(']')
        s = ''.join(str(_) for _ in s)
        return f'{s}'        