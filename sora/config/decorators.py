import functools
import warnings

next_major_version = 'v1.0'

warnings.simplefilter('always', FutureWarning)


# This decorator gets an argument that is being deprecated
def deprecated_alias(**aliases):
    def deco(f):
        @functools.wraps(f)
        def wrapper(*args, **kwargs):
            rename_kwargs(f.__name__, kwargs, aliases)
            return f(*args, **kwargs)

        return wrapper

    return deco


def rename_kwargs(func_name, kwargs, aliases):
    for alias, new in aliases.items():
        if alias in kwargs:
            if new in kwargs:
                raise TypeError('{} received both {} and {}'.format(
                    func_name, alias, new))
            warnings.warn("'{}' is deprecated and will be removed in {}; please use '{}'".
                          format(alias, next_major_version, new), FutureWarning)
            kwargs[new] = kwargs.pop(alias)


def deprecated_function(message):
    def deco(func):
        """This is a decorator which can be used to mark functions
        as deprecated. It will result in a warning being emitted
        when the function is used."""

        @functools.wraps(func)
        def new_func(*args, **kwargs):
            warnings.warn("{} is deprecated and will be removed in {}; {}".
                          format(func.__name__, next_major_version, message),
                          category=FutureWarning,
                          stacklevel=2)
            return func(*args, **kwargs)

        return new_func

    return deco
