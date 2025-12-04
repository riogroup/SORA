# sora/rings/meta.py
import numpy as np
import astropy.units as u
from sora.body.meta import PhysicalData 
from sora.config import input_tests
from astropy.coordinates import SkyCoord



__all__ = ["BaseRing"]


class BaseRing:
    """
    BaseRing — physical parameters of a planetary ring.

    Stores intrinsic (non-geometric) physical parameters such as:
    radius, width, opacity, and equivalent depth.

    Notes
    -----
    - Does not include geometric orientation (pole, inclination, etc.);
      these are handled by `RingGeometry`.
    - All attributes are stored as `PhysicalData` objects.
    """

    def __init__(self, **kwargs):
        allowed_kwargs = ['name',
                          'radius', 'radius_err',
                          'eccentricity', 'eccentricity_err',
                          'pole_orientation',
                          'normal_opacity', 'normal_opacity_err',
                          'normal_optical_depth', 'normal_optical_depth_err',
                          'radial_width', 'radial_width_err',
                          'equivalent_depth', 'equivalent_depth_err',
                          'equivalent_width', 'equivalent_width_err']
        
            
        input_tests.check_kwargs(kwargs, allowed_kwargs=allowed_kwargs)

        radius = kwargs.get("radius", None)
        if radius is None:
            self._radius = PhysicalData("Radius", None, 0.0, unit=u.km)
        else:
            self._radius = PhysicalData("Radius", radius, kwargs.get("radius_err", 0.0), unit=u.km)
        
        self._normal_opacity = PhysicalData('Normal opacity', 
                                            kwargs.get('normal_opacity', None), 
                                            kwargs.get('normal_opacity_err', 0.0))  
        
        self._normal_optical_depth = PhysicalData('Normal optical depth', 
                                                  kwargs.get('normal_optical_depth', None), 
                                                  kwargs.get('normal_optical_depth_err', 0.0), unit=u.km)
        
        self._radial_width = PhysicalData('Radial width', 
                                          kwargs.get('radial_width', None), 
                                          kwargs.get('radial_width_err', 0.0), unit=u.km)
        
        self._eccentricity = PhysicalData('Eccentricity', 
                                          kwargs.get('eccentricity', None), 
                                          kwargs.get('eccentricity_err', 0.0))
        
        self._equivalent_depth = PhysicalData('Equivalent depth', 
                                              kwargs.get('equivalent_depth', None), 
                                              kwargs.get('equivalent_depth_err', 0.0), unit=u.km)
        
        self._equivalent_width = PhysicalData('Equivalent width', 
                                              kwargs.get('equivalent_width', None), 
                                              kwargs.get('equivalent_width_err', 0.0), unit=u.km)

    @property
    def pole_orientation(self):
        """Return ring pole orientation as SkyCoord."""
        if self.geometry is None:
            return None
        return self.geometry.pole_orientation

    @pole_orientation.setter
    def pole_orientation(self, value):
        """
        Set the ring pole orientation.
        Accepts SkyCoord, or strings like "180d 30d".
        """
        if self.geometry is None:
            raise AttributeError("Ring has no geometry object associated.")

        self.geometry.pole_orientation = value

    @property
    def radius(self):
        return self._radius    
    
    @radius.setter
    def radius(self, value):
        if isinstance(value, PhysicalData):
            radius = value
            radius.name = "Radius"
        else:
            radius = PhysicalData('Radius', value)
        if radius < 0:
            raise ValueError("Radius cannot be a negative value")
        self._radius = radius        
        
    @property
    def normal_opacity(self):
        return self._normal_opacity  
      
    @normal_opacity.setter
    def normal_opacity(self, value):
        if isinstance(value, PhysicalData):
            normal_opacity = value
            normal_opacity.name = "Normal opacity"
        else:
            normal_opacity = PhysicalData('Normal opacity', value)
            
        if normal_opacity < 0 or normal_opacity > 1:
            raise ValueError("Normal opacity must be betwenn 0 and 1")
            
        self._normal_opacity = normal_opacity

    
    @property
    def normal_optical_depth(self):
        return self._normal_optical_depth    
    
    @normal_optical_depth.setter
    def normal_optical_depth(self, value):
        if isinstance(value, PhysicalData):
            normal_optical_depth = value
            normal_optical_depth.name = "Normal optical depth"
        else:
            normal_optical_depth = PhysicalData('Normal optical depth', value)
            
        if normal_optical_depth < 0:
            raise ValueError("Normal optical depth must be positive")
            
        self._normal_optical_depth = normal_optical_depth
        
    @property
    def radial_width(self):
        return self._radial_width    
    
    @radial_width.setter
    def radial_width(self, value):
        if isinstance(value, PhysicalData):
            radial_width = value
            radial_width.name = "Radial width"
        else:
            radial_width = PhysicalData('Radial width', value, unit=u.km)
            
        if radial_width < 0:
            raise ValueError("Radial width must be positive")
            
        self._radial_width = radial_width
        
    @property
    def eccentricity(self):
        return self._eccentricity    
    
    @eccentricity.setter
    def eccentricity(self, value):
        if isinstance(value, PhysicalData):
            eccentricity = value
            eccentricity.name = "Eccentricity"
        else:
            eccentricity = PhysicalData('Eccentricity', value)
            
        if eccentricity < 0 or eccentricity > 1:
            raise ValueError("Eccentricity must be betwenn 0 and 1")
            
        self._eccentricity = eccentricity

    @property
    def equivalent_depth(self):
        return self._equivalent_depth   
     
    @equivalent_depth.setter
    def equivalent_depth(self, value):
        if isinstance(value, PhysicalData):
            equivalent_depth = value
            equivalent_depth.name = "Equivalent depth"
        else:
            equivalent_depth = PhysicalData('Equivalent depth', value)
            
        if equivalent_depth < 0:
            raise ValueError("Equivalent depth must be positive")
            
        self._equivalent_depth = equivalent_depth

    @property
    def equivalent_width(self):
        return self._equivalent_width 
    
    @equivalent_width.setter
    def equivalent_width(self, value):
        if isinstance(value, PhysicalData):
            equivalent_width = value
            equivalent_width.name = "Equivalent width"
        else:
            equivalent_width = PhysicalData('Equivalent width', value)
            
        if equivalent_width < 0:
            raise ValueError("Equivalent width must be positive")
            
        self._equivalent_width = equivalent_width

    def summary(self):
        """Return string summary of all defined parameters."""
        fields = ["radius", "radial_width", "normal_opacity", "normal_optical_depth",
                  "equivalent_depth", "equivalent_width", "eccentricity"]
        out = ["Ring physical parameters:"]
        for f in fields:
            val = getattr(self, f, None)
            if val is not None:
                out.append(f"  - {f.replace('_',' ').title()}: {val}")
        return "\n".join(out)

    def __str__(self):
        return self.summary()
