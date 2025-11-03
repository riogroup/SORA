from .core import *
from .utils import *
from .core import LightCurve
from .model import attach_to_lightcurve_class as _attach_models
_attach_models(LightCurve)

__all__ = ['LightCurve']
