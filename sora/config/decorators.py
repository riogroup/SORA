import functools
import warnings

next_major_version = 'v1.0'

warnings.simplefilter('always', FutureWarning)


def deprecated_alias(**aliases):
    """Create a decorator that renames deprecated keyword arguments.

    Parameters
    ----------
    **aliases
        Mapping from deprecated keyword names to their replacement names.

    Returns
    -------
    decorator : `function`
        Decorator that updates deprecated keyword arguments and emits a
        `FutureWarning` when an alias is used.
    """
    def deco(f):
        @functools.wraps(f)
        def wrapper(*args, **kwargs):
            rename_kwargs(f.__name__, kwargs, aliases)
            return f(*args, **kwargs)

        return wrapper

    return deco


def rename_kwargs(func_name, kwargs, aliases):
    """Rename deprecated keyword arguments in place.

    Parameters
    ----------
    func_name : `str`
        Name of the function that received the keyword arguments.
    kwargs : `dict`
        Keyword arguments to inspect and update.
    aliases : `dict`
        Mapping from deprecated keyword names to their replacement names.

    Raises
    ------
    TypeError
        If both a deprecated keyword and its replacement are provided.
    """
    for alias, new in aliases.items():
        if alias in kwargs:
            if new in kwargs:
                raise TypeError('{} received both {} and {}'.format(
                    func_name, alias, new))
            warnings.warn("'{}' is deprecated and will be removed in {}; please use '{}'".
                          format(alias, next_major_version, new), FutureWarning)
            kwargs[new] = kwargs.pop(alias)


def deprecated_function(message):
    """Create a decorator that marks a function as deprecated.

    Parameters
    ----------
    message : `str`
        Additional deprecation message shown after the function name and
        removal version.

    Returns
    -------
    decorator : `function`
        Decorator that emits a `FutureWarning` whenever the wrapped function is
        called.
    """
    def deco(func):
        @functools.wraps(func)
        def new_func(*args, **kwargs):
            warnings.warn("{} is deprecated and will be removed in {}; {}".
                          format(func.__name__, next_major_version, message),
                          category=FutureWarning,
                          stacklevel=2)
            return func(*args, **kwargs)

        return new_func

    return deco
