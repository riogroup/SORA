def check_kwargs(input_kwargs, allowed_kwargs, raise_error=True):
    """Test whether input keyword names are allowed.

    Parameters
    ----------
    input_kwargs : `dict`, `list`
        Dictionary or list with the input keyword names.

    allowed_kwargs : `list`
        List with the allowed keys.

    raise_error : `bool`
        If `True`, raise an error when unexpected keyword names are found. If
        `False`, return the unexpected keyword names instead.

    Returns
    -------
    not_allowed : `list`
        List with unexpected keyword names. Returned only when
        ``raise_error=False``.

    Raises
    ------
    TypeError
        If ``raise_error=True`` and at least one keyword name is not allowed.
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
    """Test and convert an attribute to an expected type.

    Parameters
    ----------
    attr
        Attribute value to test.
    typ : `type`
        Expected type or callable used to convert `attr`.
    name : `str`
        Name of the variable used in the error message.

    Returns
    -------
    attr
        `attr` converted by `typ`.

    Raises
    ------
    TypeError
        If `attr` cannot be converted by `typ`.
    """
    try:
        return typ(attr)
    except:
        raise TypeError('"{}" keyword must be a {} object'.format(name, typ))


class SelectDefault:
    """Validate parameters against explicit values and named defaults."""

    def __init__(self, instance, defaults: dict):
        """Initialize the selector.

        This class is not meant to be used directly by users. It facilitates
        testing parameters that have default values. It accepts either a string
        key for a default parameter or a parameter with an allowed type.

        Parameters
        ----------
        instance : `type`
            Allowed type for explicit parameter values.

        defaults : `dict`
            Dictionary with the keys and values for default parameters.

        Raises
        ------
        TypeError
            If a default key is not a string or a default value does not have
            the allowed type.

        """
        for key, value in defaults.items():
            if not isinstance(key, str):
                raise TypeError("key '{}' must be a string".format(key))
            if not isinstance(value, instance):
                raise TypeError("{} is not an allowed type {}".format(key, instance))
        self.instance = instance
        self.allowed_keys = defaults

    def get_default(self, value):
        """Return a default or validate an explicit parameter value.

        Parameters
        ----------
        value
            Default key or explicit parameter value.

        Returns
        -------
        value
            The selected default value or the validated explicit value.

        Raises
        ------
        ValueError
            If `value` is neither a known default key nor an allowed explicit
            value.
        """
        if isinstance(value, str):
            value = self.allowed_keys.get(value, value)
        if not isinstance(value, self.instance):
            raise ValueError(
                "{} is not an allowed parameter '{}' or a type '{}'".format(value, self.allowed_keys.keys(), self.instance))
        return value
