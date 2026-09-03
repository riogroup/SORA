from importlib import import_module
from typing import TYPE_CHECKING


__all__ = [
    "Body",
    "EphemHorizons",
    "EphemJPL",
    "EphemKernel",
    "EphemPlanete",
    "LightCurve",
    "Observer",
    "Occultation",
    "Spacecraft",
    "Star",
    "__version__",
]

_LAZY_IMPORTS = {
    "Body": (".body", "Body"),
    "EphemHorizons": (".ephem", "EphemHorizons"),
    "EphemJPL": (".ephem", "EphemJPL"),
    "EphemKernel": (".ephem", "EphemKernel"),
    "EphemPlanete": (".ephem", "EphemPlanete"),
    "LightCurve": (".lightcurve", "LightCurve"),
    "Observer": (".observer", "Observer"),
    "Occultation": (".occultation", "Occultation"),
    "Spacecraft": (".observer", "Spacecraft"),
    "Star": (".star", "Star"),
    "__version__": (".version", "version"),
}


if TYPE_CHECKING:
    from .body import Body
    from .ephem import EphemHorizons, EphemJPL, EphemKernel, EphemPlanete
    from .lightcurve import LightCurve
    from .observer import Observer, Spacecraft
    from .occultation import Occultation
    from .star import Star
    from .version import version as __version__


def __getattr__(name):
    """Load public objects only when they are first requested."""
    try:
        module_name, attribute_name = _LAZY_IMPORTS[name]
    except KeyError:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}") from None

    value = getattr(import_module(module_name, __name__), attribute_name)
    globals()[name] = value
    return value


def __dir__():
    return sorted(set(globals()) | set(__all__))
