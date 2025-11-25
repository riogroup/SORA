# sora/rings/__init__.py
"""
sora.rings
==========

Ring module — handles physical, geometric, and observational properties
of planetary rings.

Public classes and functions:
-----------------------------
- Ring : main class (physical + geometric data)
- BaseRing : low-level physical data container
- RingGeometry : orientation and projection handler
- compute_local_properties : derive ring parameters from occultation chords
"""

from .core import Ring
from .meta import BaseRing
from .geometry import RingGeometry
from .properties import compute_local_properties
from .utils import *

__all__ = [
    "Ring",
    "BaseRing",
    "RingGeometry",
    "compute_local_properties"
]
