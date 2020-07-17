def check_kwargs(input_kwargs, allowed_kwargs, raise_error=True):
    """ Tests if the input kwargs are allowed
    
    Parameters:
        input_kwargs: Dictionary or list with the input values
        allowed_kwargs: list with the allowed keys
        raise_error: if True, it will raises an error,
            if False, it returns the keys not allowed.
    """
    not_allowed = [i for i in input_kwargs if i not in allowed_kwargs]
    if raise_error:
        if len(not_allowed) > 0:
            allowed_kwargs.sort()
            raise TypeError("function got an unexpected keyword argument {}\n"
                            "Available kwargs are: {}".format(not_allowed, allowed_kwargs))
    else:
        return not_allowed

