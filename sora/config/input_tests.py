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


class SelectDefault:
    def __init__(self, instance, defaults: dict):
        """ This class is not meant to be used by the user.

        The goal of this class is to facilitate testing parameters
        that have default values. So when using, the class will test if
        the user passed a string with the key for a default parameter,
        or a parameter with an allowed type.

        Parameters
        ----------
        instance : type
            The instance of the allowed parameters

        defaults: dict
            A dictionary with the keys and values for default parameters.

        """
        for key, value in defaults.items():
            if not isinstance(key, str):
                raise TypeError("key '{}' must be a string".format(key))
            if not isinstance(value, instance):
                raise TypeError("{} is not an allowed type {}".format(key, instance))
        self.instance = instance
        self.allowed_keys = defaults

    def get_default(self, value):
        if isinstance(value, str):
            value = self.allowed_keys.get(value, value)
        if not isinstance(value, self.instance):
            raise ValueError(
                "{} is not an allowed parameter '{}' or a type '{}'".format(value, self.allowed_keys.keys(), self.instance))
        return value
