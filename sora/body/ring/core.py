from sora.config import input_tests
from sora.body.meta import PhysicalData
from sora.body import Body
from sora.ephem.meta import BaseEphem
from astropy.coordinates import SkyCoord
from astropy.time import Time
from .meta import BaseRing
import astropy.units as u
import numpy as np
import warnings
from sora.body.ring.utils import calc_coef_projecao, project_to_ring_plane

class Ring(BaseRing):    
    def __init__(self, ephem, **kwargs):
        """ Class that contains and manages the information of a ring.

        Parameters:
            ring_id : `str`, required 
                Identifier of the ring. 

            radius :  `float`, `int`, `astropy.quantity.Quantity`, required
                Mean radial distance of the ring from the central body, in km.

            radius_err : `float`, `int`
                Uncertainty associated with the ring radius.

            eccentricity : `float`
                Orbital eccentricity of the ring.

            eccentricity_err : `float`
                Uncertainty associated with the eccentricity.

            pole_orientation : `str`, `astropy.coordinates.SkyCoord`
                The Pole coordinates of the ring. It can be a `SkyCoord object` or a
                string in the format ``'hh mm ss.ss +dd mm ss.ss'``.
               
            normal_opacity : `float`, `int`
                Normal opacity of the ring, the opacity corrected for the ring opening angle.
                           
            normal_opacity_err : `float`, `int`
                Uncertainty associated with the normal opacity.

            normal_optical_depth : `float`, `int`
                Normal optical depth of the ring, the optical depth corrected for the ring opening angle.

            normal_optical_depth_err : `float`, `int`
                Uncertainty associated with the normal optical depth.

            radial_width : `float`, `int`, `astropy.quantity.Quantity`
                Radial width of the ring, in kilometers.

            radial_width_err (float): 
                Uncertainty in the measured radial width, in kilometers.

            equivalent_depth : `float`, `int`, `astropy.quantity.Quantity` 
                Equivalent depth of the ring, defined as the integral of the normal optical depth over the ring radial width, in kilometers.

            equivalent_depth_err : `float`, `int`, `astropy.quantity.Quantity` 
                Uncertainty in the equivalent depth, in kilometers.

            equivalent_width : `float`, `int`, `astropy.quantity.Quantity` 
                Equivalent width of the ring, defined as the integral of the normal optical depth across the radial profile, in kilometers.

            equivalent_width_err : `float`, `int`, `astropy.quantity.Quantity` 
                Uncertainty in the equivalent width, in kilometers.

        """
        
        allowed_kwargs = ['ring_id',
                          'radius', 'radius_err',
                          'eccentricity', 'eccentricity_err',
                          'pole_orientation',
                          'normal_opacity', 'normal_opacity_err',
                          'normal_optical_depth', 'normal_optical_depth_err',
                          'radial_width', 'radial_width_err',
                          'equivalent_depth', 'equivalent_depth_err',
                          'equivalent_width', 'equivalent_width_err']
        
            
        input_tests.check_kwargs(kwargs, allowed_kwargs=allowed_kwargs)

        self.ephem = ephem        
        self.ring_id = kwargs.get('ring_id')
        
        if not 'pole_orientation' in kwargs:
            raise ValueError(f'Pole orientation is not defined for ring {self.ring_id}') 
        else:
            pole = SkyCoord(kwargs.get('pole_orientation'))        
            self._pole_orientation = pole.icrs

        if not 'radius' in kwargs:
            raise ValueError(f'Radius is not defined for ring {self.ring_id}')
        else:            
            self._radius = PhysicalData('Radius', 
                                        kwargs.get('radius'), 
                                        kwargs.get('radius_err', 0.0), unit=u.km)    
        
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
        

    def get_ring_orientation(self, time, observer='geocenter'):
        """
        Computes the ring's position and opening angles as seen by a given observer at a specified time.

        Parameters
        ----------
        time : `str`, `astropy.time.Time`
            Reference time to calculate the position. It can be a string
            in the ISO format (yyyy-mm-dd hh:mm:ss.s) or an astropy Time object.

        observer : `str`, `sora.Observer`, `sora.Spacecraft`
            IAU code of the observer (must be present in given list of kernels),
            a SORA observer object or a string: ['geocenter', 'barycenter']

        Returns:
            tuple : astropy.units.Quantity, astropy.units.Quantity
                Ring position angle and opening angle, both in degrees.
        """

        time = Time(time)
        pos = self.ephem.get_position(time, observer=observer)
        ring_position_angle = pos.position_angle(self._pole_orientation).rad*u.rad
        ring_opening_angle = np.arcsin(-(np.sin(self._pole_orientation.dec)*np.sin(pos.dec) + 
                                          np.cos(self._pole_orientation.dec)*np.cos(pos.dec)*np.cos(self._pole_orientation.ra-pos.ra)
                                         )
                                       )
        return ring_position_angle.to('deg'), ring_opening_angle.to('deg')
    
    def to_ring_plane(self, f, g, time, center_f=0, center_g=0):
        """
        Convert sky-plane coordinates (f, g) to ring-plane (x, y),
        computing internally the projection coefficients.

        Parameters
        ----------
        f, g : array_like
            Sky-plane coordinates [km]
        time : astropy.time.Time
            Time of observation.
        ephem : Ephem object
            Ephemeris of the body.
        center_f, center_g : float
            Sky-plane offset [km]

        
        Returns
        -------
        x, y : ndarray
            Coordinates in the ring plane [km]
        """
        pos = self.ephem.get_position(time)
        P, B = self.get_ring_orientation(time)
        earth_pole = SkyCoord('12h00m00s +90d00m00s')

        coef, coef_polo = calc_coef_projecao(pos, self.pole_orientation, B, P, earth_pole)
        x, y = project_to_ring_plane(f, g, coef, coef_polo, ksi_0=center_f, eta_0=center_g)

        return x, y
    
    def __str__(self):

        """ String representation of the Ring class
        """

        out = []
        out.append('Pole orientation (ICRS):\n    RA: {:.2f}\n    DEC: {:.2f}\n'.format(self._pole_orientation.icrs.ra,
                                                                                        self._pole_orientation.icrs.dec)
                  )
        #out.append('Pole orientation (Barycentric Mean Ecliptic):\n    lat: {:.4f}\n    lon: {:.4f}\n'.format(self.pole_orientation.barycentricmeanecliptic.lat,
        #                                                                                  self.pole_orientation.barycentricmeanecliptic.lon)
        #          )
        
        out.append(self._radius.__str__())        
        out.append(self._normal_opacity.__str__())
        out.append(self._normal_optical_depth.__str__())
        out.append(self._radial_width.__str__())
        out.append(self._eccentricity.__str__())
        out.append(self._equivalent_width.__str__())  
        out.append(self._equivalent_depth.__str__())
        
        return ''.join(out)
