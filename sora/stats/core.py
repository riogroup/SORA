import numpy as np


class Parameter:
    '''A collection of values in an object to be used in a fitting procedure.
       The structure of this class heavely borrows its structure from parameters.py (lmfit)
    '''

    def __init__(self, name, value=None, minval=-np.inf, maxval=np.inf, free=True, std=np.inf, initial_value=None):
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
                print(parobj)
                print(parobj.name)
            
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


    # THIS METHOD IS TEMPORARY AND NEEDS A BETTER SOLUTION
    def summary(self):
        '''
        Prints a summary of Parameters object.

        Returns
        -------
        None

        '''
        
        keys = []
        for key in self.valuesdict().keys():
            keys.append(key)

        labels = [' Name     ', ' Value       ', ' Lower Bound ', ' Upper Bound ', '    Free    ']
        keylenghts = [ len(k) for k in keys ]
        if max(keylenghts) > len(labels[0]):

            textmasklabels = '|{:'+str(max(keylenghts)+2)+'}|{:13}|{:13}|{:13}|{:13}|'
            dash = '-'*13
            dashmask = '-'*(max(keylenghts)+2)
            print('+{}+{}+{}+{}+{}+'.format(dashmask,dash,dash,dash,dash))
            print(textmasklabels.format(labels[0], labels[1], labels[2], labels[3], labels[4]))
            print('+{}+{}+{}+{}+{}+'.format(dashmask,dash,dash,dash,dash))
            textmask = '|{:'+str(max(keylenghts)+2)+'}|{:13f}|{:13f}|{:13f}|    {:9}|'
            for key in keys:
                print(textmask.format(self[key].name, self[key].value, self[key].minval, self[key].maxval, 'True' if self[key].free else 'False'))
            print('+{}+{}+{}+{}+{}+'.format(dashmask,dash,dash,dash,dash))
        else:
            dash = '-'*13
            print('+{}+{}+{}+{}+{}+'.format(dash,dash,dash,dash,dash))
            print('|{:13}|{:13}|{:13}|{:13}|{:13}|'.format(labels[0], labels[1], labels[2], labels[3], labels[4]))
            print('+{}+{}+{}+{}+{}+'.format(dash,dash,dash,dash,dash))
            for key in keys:
                print('|{:13}|{:13f}|{:13f}|{:13f}|    {:9}|'.format(self[key].name, self[key].value, self[key].minval, self[key].maxval, 'True' if self[key].free else 'False'))
            print('+{}+{}+{}+{}+{}+'.format(dash,dash,dash,dash,dash))

            