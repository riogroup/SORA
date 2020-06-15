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
