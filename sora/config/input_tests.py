def check_kwargs(input_kwargs, allowed_kwargs, raise_error=True):
    """Tests if the input `**kwargs` are allowed.

    Parameters
    ----------
    input_kwargs : `dict`, `list`
        Dictionary or list with the input values.

    allowed_kwargs : `list`
        List with the allowed keys.

    raise_error : `bool`
        Raises error if ``True``. If ``False``, it will raise the error listing
        the not allowed keys.
    """
    not_allowed = [i for i in input_kwargs if i not in allowed_kwargs]
    if raise_error:
        if len(not_allowed) > 0:
            allowed_kwargs.sort()
            raise TypeError("function got an unexpected keyword argument {}\n"
                            "Available kwargs are: {}".format(not_allowed, allowed_kwargs))
    else:
        return not_allowed


def test_attr(attr, typ, name):
    """Tests attribute.

    This function tests if the attribute ``attr`` belongs to the type ``typ``,
    if not it gives an error message informing the name of the variable.
    """
    try:
        return typ(attr)
    except:
        raise TypeError('"{}" keyword must be a {} object'.format(name, typ))
