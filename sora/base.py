import numpy as np


class Parameter:
    '''A collection of values in an object to be used in a fitting procedure.
       The structure of this class heavely borrows its structure from parameters.py (lmfit)
    '''

    def __init__(self, name, value=None, std=inf, initial_value=None, minval=-inf, maxval=inf, free=True):
        '''
        Constructor method for the Parameter object.

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

        if isclose( self.minval, self.maxval, atol=1e-15, rtol=1e-15):
            raise ValueError(f'Parameter {self.name} has `minval` equal to maxval')
        
        if self.minval <= value <= self.maxval:
            self.value = value
        else:
            raise ValueError(f'Value `{value}` at parameter `{self.name}` is outside the provided bounds `({self.minval}, {self.maxval})`')
            
            
        
    def update(self, **kwargs):
        for key, value in kwargs.items():
            setattr(self, key, value)