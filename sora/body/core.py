from sora.config import input_tests
from .utils import search_sbdb, apparent_magnitude
from .meta import BaseBody, PhysicalData
import astropy.units as u
from astropy.coordinates import SkyCoord, Longitude, Latitude, get_sun
from astropy.time import Time
from astroquery.jplhorizons import Horizons
import numpy as np
import warnings


__all__ = ['Body']


class Body(BaseBody):

    def __init__(self, name, mode='sbdb', orbit_class=None, spkid=None, **kwargs):
        """ Class that contains and manage the information of the body

        Parameters:
            name (str): The name of the object. It can be the used spkid or designation number
                to query the SBDB. In this case, the name is case insensitive. (Required)
            mode (str): If mode='sbdb', a query on the Small-Body DataBase is made
                to identify physical parameter of the object. If mode='local', it uses the
                user parameters. In this case, "orbit_class" and "spkid" parameters may be required.
            ephem (EphemKernel, EphemHorizons, EphemJPL, EphemPlanete): An Ephem Class that
                contains information about the ephemeris. It can be "horizons" to automatically
                defined an EphemHorizons object or a list of kernels to automatically define an
                EphemKernel object.

        Parameters that are returned from the Small-Body DataBase if mode='sbdbd'. These are
            the physical paramaters the user can give to the object. If a query is made and
            user gives a parameter, the parameter given by the user is defined in the Body object:

            orbit_class (str): It defines the Orbital class of the body. It can be 'TNO',
                'Satellite', 'Centaur', 'comet', 'asteroid', 'trojan', 'neo' and 'planet'.
                It is important for a better characterization of the object.
                It is required if mode='local', else it will be replaced by the queried value.
            spkid (str, number): If mode='local', the user must give a spkid or an ephem
                which has the spkid parameter.
            albedo (number): The albedo of the object.
            H (number): The absolute magnitude.
            G (number): The phase slope.
            diameter (number, Quantity): The dimater of the object, in km.
            density: (number, Quantity): The density of the object, in g/cm^3.
            GM (number, Quantity): The Standard Gravitational Parameter, in km^3/s^2.
            rotation (number, Quantity): The Rotation of the object, in hours.
            pole (str, SkyCoord): The Pole coordinates of the object. It can be a SkyCoord
                object or a string in the format 'hh mm ss.ss +dd mm ss.ss'
            BV (number): The B-V color.
            UB (number): The U-B color.
            smass (str): The spectral type in SMASS classification.
            tholen (str): The spectral type in Tholen classification.
        """
        allowed_kwargs = ["albedo", "H", "G", "diameter", "density", "GM", "rotation", "pole", "BV", "UB", "smass", "tholen",
                          "ephem"]
        input_tests.check_kwargs(kwargs, allowed_kwargs=allowed_kwargs)
        self._shared_with = {'ephem': {}, 'occultation': {}}
        self.orbit_class = orbit_class
        if spkid:
            self.spkid = spkid
        modes = {'sbdb': self.__from_sbdb, 'local': self.__from_local}
        if mode not in modes:
            raise ValueError("'mode' parameter can only be 'sbdb' or 'local'")
        modes[mode](name=name)

        # set the physical parameters based on the kwarg name.
        if 'smass' in kwargs:
            self.spectral_type['SMASS']['value'] = kwargs.pop('smass')
        if 'tholen' in kwargs:
            self.spectral_type['Tholen']['value'] = kwargs.pop('tholen')
        if spkid:
            self.spkid = spkid
        for key in kwargs:
            setattr(self, key, kwargs[key])
        self._shared_with['ephem']['search_name'] = self._search_name
        self._shared_with['ephem']['id_type'] = self._id_type

    def __from_sbdb(self, name):
        """Search the object in SBDB and define its physical parameters

        Parameters:
            name (str): The name, spkid or designation number of the Small Body.
        """
        sbdb = search_sbdb(name)
        self.meta_sbdb = sbdb
        self.name = sbdb['object']['fullname']
        self.shortname = sbdb['object'].get('shortname', self.name)
        self.orbit_class = sbdb['object']['orbit_class']['name']

        pp = sbdb['phys_par']  # get the physical parameters (pp) of the sbdb

        self.albedo = PhysicalData('Albedo', pp.get('albedo'), pp.get('albedo_sig'), pp.get('albedo_ref'), pp.get('albedo_note'))
        self.H = PhysicalData('Absolute Magnitude', pp.get('H'), pp.get('H_sig'), pp.get('H_ref'), pp.get('H_note'), unit=u.mag)
        self.G = PhysicalData('Phase Slope', pp.get('G'), pp.get('G_sig'), pp.get('G_ref'), pp.get('G_note'))
        self.diameter = PhysicalData('Diameter', pp.get('diameter'), pp.get('diameter_sig'), pp.get('diameter_ref'),
                                     pp.get('diameter_note'), unit=u.km)
        self.density = PhysicalData('Density', pp.get('density'), pp.get('density_sig'), pp.get('density_ref'),
                                    pp.get('density_note'), unit=u.g/u.cm**3)
        self.GM = PhysicalData('Standard Gravitational Parameter', pp.get('GM'), pp.get('GM_sig'), pp.get('GM_ref'),
                               pp.get('GM_note'), unit=u.km**3/u.s**2)
        self.rotation = PhysicalData('Rotation', pp.get('rot_per'), pp.get('rot_per_sig'), pp.get('rot_per_ref'),
                                     pp.get('rot_per_note'), unit=u.h)
        self.BV = PhysicalData('B-V color', pp.get('BV'), pp.get('BV_sig'), pp.get('BV_ref'), pp.get('BV_note'))
        self.UB = PhysicalData('U-B color', pp.get('UB'), pp.get('UB_sig'), pp.get('UB_ref'), pp.get('UB_note'))
        if 'pole' in pp:
            self.pole = SkyCoord(pp['pole'].replace('/', ' '), unit=('deg', 'deg'))
            self.pole.ra.uncertainty = Longitude(pp['pole_sig'].split('/')[0], unit=u.deg)
            self.pole.dec.uncertainty = Latitude(pp['pole_sig'].split('/')[1], unit=u.deg)
            self.pole.reference = pp['pole_ref'] or ""
            self.pole.notes = pp['pole_note'] or ""
        else:
            self.pole = None
        self.spectral_type = {
            "SMASS": {"value": pp.get('spec_B'), "reference": pp.get('spec_B_ref'), "notes": pp.get('spec_B_note')},
            "Tholen": {"value": pp.get('spec_T'), "reference": pp.get('spec_T_ref'), "notes": pp.get('spec_T_note')}}
        self.spkid = sbdb['object']['spkid']
        self._des_name = sbdb['object']['des']
        self.discovery = "Discovered {} by {} at {}".format(sbdb['discovery'].get('date'), sbdb['discovery'].get('who'),
                                                            sbdb['discovery'].get('location'))

    def __from_local(self, name):
        """Define Body object with default values for mode="local"

        Parameters:

        """
        self.name = name
        self.shortname = name
        orbit_classes = {'tno': 'TransNeptunian Object', 'satellite': 'Natural Satellite', 'centaur': 'Centaur',
                         'comet': 'Comet', 'asteroid': 'Main-belt Asteroid', 'trojan': 'Jupiter Trojan',
                         'neo': 'Near-Earth Object', 'planet': "Planet"}
        self.orbit_class = orbit_classes.get(self.orbit_class.lower())
        if not self.orbit_class:
            raise ValueError(("'orbit_class' must be one of the following: {}".
                             format([k.capitalize() for k in orbit_classes.keys()])))
        self.albedo = None
        self.H = None
        self.G = None
        self.diameter = None
        self.density = None
        self.GM = None
        self.rotation = None
        self.pole = None
        self.BV = None
        self.UB = None
        self.spectral_type = {"SMASS": {"value": None, "reference": None, "notes": None},
                              "Tholen": {"value": None, "reference": None, "notes": None}}
        self.discovery = ""

    def get_pole_position_angle(self, time):
        """ Returns the pole position angle and aperture angle relative to the geocenter

        Parameters:
            time (str, Time): Time from which to calculate the position.

        Returns:
            position_angle (float): Position angle of the object pole, in degrees
            aperture_angle (float): Apeture angle of the object pole, in degrees
        """
        time = Time(time)
        pole = self.pole
        if np.isnan(pole.ra):
            raise ValueError("Pole coordinates are not defined")
        obj = self.ephem.get_position(time)
        position_angle = obj.position_angle(pole).rad*u.rad
        aperture_angle = np.arcsin(
            -(np.sin(pole.dec)*np.sin(obj.dec) +
              np.cos(pole.dec)*np.cos(obj.dec)*np.cos(pole.ra-obj.ra))
            )
        return position_angle.to('deg'), aperture_angle.to('deg')

    def apparent_magnitude(self, time):
        """ Calculates the Object Apparent Magnitude

        Parameters:
            time (str, Time): Reference time to calculate the object aparent magnitude.

        Returns:
            ap_mag (float): Object apparent magnitude
        """
        time = Time(time)

        if np.isnan(self.H) or np.isnan(self.G):
            warnings.warn('H and/or G is not defined for {}. Searching into JPL Horizons service'.format(self.shortname))
            obj = Horizons(id=self._search_name, id_type=self._id_type, location='geo', epochs=time.jd)
            eph = obj.ephemerides(extra_precision=True)
            if 'H' in eph.keys():
                self.H = eph['H'][0]
                self.H.reference = "JPL Horizons"
                self.G = eph['G'][0]
                self.G.reference = "JPL Horizons"
            if len(eph['V']) == 1:
                return eph['V'][0]
            else:
                return eph['V'].tolist()

        else:
            obs_obj = self.ephem.get_position(time)
            obs_sun = get_sun(time)
            sun_obj = SkyCoord(obs_obj.cartesian - obs_sun.cartesian)
            sun_obj.representation_type = 'spherical'

            # Calculates the phase angle between the 2-vectors
            unit_vector_1 = -obs_obj.cartesian.xyz / np.linalg.norm(obs_obj.cartesian.xyz)
            unit_vector_2 = -sun_obj.cartesian.xyz / np.linalg.norm(sun_obj.cartesian.xyz)
            dot_product = np.dot(unit_vector_1, unit_vector_2)
            phase = np.arccos(dot_product).to(u.deg).value
            return apparent_magnitude(self.H.value, self.G.value, obs_obj.distance.to(u.AU).value,
                                      sun_obj.distance.to(u.AU).value, phase)

    def __str__(self):
        from .values import smass, tholen
        out = []
        out.append('#'*79 + '\n{:^79s}\n'.format(self.name) + '#'*79 + '\n')
        out.append('Object Orbital Class: {}\n'.format(self.orbit_class))
        if self.spectral_type['Tholen']['value'] or self.spectral_type['SMASS']['value']:
            out += 'Spectral Type:\n'
            value = self.spectral_type['SMASS']['value']
            if value:
                out.append('    SMASS: {}  [Reference: {}]\n'.format(value, self.spectral_type['SMASS']['reference']))
            value = self.spectral_type['Tholen']['value']
            if value:
                out.append('    Tholen: {} [Reference: {}]\n'.format(value, self.spectral_type['Tholen']['reference']))
            out += " "*7 + (smass.get(self.spectral_type['SMASS']['value']) or
                            tholen.get(self.spectral_type['Tholen']['value'])) + "\n"
        out.append(self.discovery)

        out.append('\n\nPhysical parameters:\n')
        out.append(self.diameter.__str__())
        out.append(self.mass.__str__())
        out.append(self.density.__str__())
        out.append(self.rotation.__str__())
        if not np.isnan(self.pole.ra):
            out.append('Pole\n    RA:{} +/- {}\n    DEC:{} +/- {}\n    Reference: {}, {}\n'.format(
                       self.pole.ra.__str__(), self.pole.ra.uncertainty.__str__(), self.pole.dec.__str__(),
                       self.pole.dec.uncertainty.__str__(), self.pole.reference, self.pole.notes))
        out.append(self.H.__str__())
        out.append(self.G.__str__())
        out.append(self.albedo.__str__())
        out.append(self.BV.__str__())
        out.append(self.UB.__str__())
        if hasattr(self, 'ephem'):
            out.append('\n' + self.ephem.__str__() + '\n')
        return ''.join(out)
