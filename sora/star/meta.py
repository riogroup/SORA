import astropy.units as u
import numpy as np
from astropy.time import Time


class MetaStar:
    """Mixin with astrometric properties used by `sora.star.Star`."""

    @property
    def ra(self):
        """Returns the right ascension of the star."""
        return self._attributes['RA']

    @ra.setter
    def ra(self, value):
        """Defines the right ascension.

        Parameters
        ----------
        value : `int`, `float`, `astropy.units.Quantity`
            Right ascension. Numeric values are interpreted as hour angle.
        """
        from astropy.coordinates import Longitude
        self._attributes['RA'] = Longitude(value, unit=u.hourangle)

    @property
    def dec(self):
        """Returns the declination of the star."""
        return self._attributes['DEC']

    @dec.setter
    def dec(self, value):
        """Defines the declination.

        Parameters
        ----------
        value : `int`, `float`, `astropy.units.Quantity`
            Declination, in degrees when unitless.
        """
        from astropy.coordinates import Latitude
        self._attributes['DEC'] = Latitude(value, unit=u.deg)

    @property
    def parallax(self):
        """Returns the parallax of the star."""
        if self.bjones:
            return self._attributes['bjones_par']
        else:
            return self._attributes.get('PAR', 0*u.mas)

    @parallax.setter
    def parallax(self, value):
        """Defines the parallax.

        Parameters
        ----------
        value : `int`, `float`, `astropy.units.Quantity`
            Parallax, in mas when unitless.
        """
        par = u.Quantity(value, unit=u.mas)
        if par <= 0*u.mas:
            par = 0*u.mas
        self._attributes['PAR'] = par

    @property
    def distance(self):
        """Returns the stellar distance inferred from parallax."""
        from astropy.coordinates import Distance
        if self.parallax > 0*u.mas:
            return Distance(parallax=self.parallax, allow_negative=False)
        else:
            raise ValueError('SORA is not able to determine distance from paralax {}'.format(self.parallax))

    @property
    def coord(self):
        """Returns the ICRS coordinate of the star."""
        from astropy.coordinates import SkyCoord
        try:
            return SkyCoord(self.ra, self.dec, self.distance)
        except ValueError:
            return SkyCoord(self.ra, self.dec)

    @property
    def pmra(self):
        """Returns the proper motion in right ascension of the star."""
        return self._attributes.get('PMRA', 0*u.mas/u.year)

    @pmra.setter
    def pmra(self, value):
        """Defines the proper motion in RA*cos(DEC).

        Parameters
        ----------
        value : `int`, `float`, `astropy.units.Quantity`
            Proper motion in RA*cos(DEC), in mas/year when unitless.
        """
        self._attributes['PMRA'] = u.Quantity(value, unit=u.mas/u.year)

    @property
    def pmdec(self):
        """Returns the proper motion in declination of the star."""
        return self._attributes.get('PMDEC', 0*u.mas/u.year)

    @pmdec.setter
    def pmdec(self, value):
        """Defines the proper motion in declination.

        Parameters
        ----------
        value : `int`, `float`, `astropy.units.Quantity`
            Proper motion in declination, in mas/year when unitless.
        """
        self._attributes['PMDEC'] = u.Quantity(value, unit=u.mas/u.year)

    @property
    def rad_vel(self):
        """Returns the radial velocity of the star."""
        return self._attributes.get('RAD_VEL', 0*u.km/u.s)

    @rad_vel.setter
    def rad_vel(self, value):
        """Defines the radial velocity.

        Parameters
        ----------
        value : `int`, `float`, `astropy.units.Quantity`
            Radial velocity, in km/s when unitless.
        """
        self._attributes['RAD_VEL'] = u.Quantity(value, unit=u.km/u.s)

    @property
    def epoch(self):
        """Returns the epoch of the star position."""
        return self._attributes['EPOCH']

    @epoch.setter
    def epoch(self, value):
        """Defines the epoch of the star position.

        Parameters
        ----------
        value : `str`, `int`, `float`, `astropy.time.Time`
            Epoch value accepted by `astropy.time.Time`.
        """
        self._attributes['EPOCH'] = Time(value)

    @property
    def bjones(self):
        """Returns whether the Bailer-Jones distance is used."""
        return self._bjones

    @bjones.setter
    def bjones(self, value):
        """Enables or disables use of Bailer-Jones distance estimates.

        Parameters
        ----------
        value : `bool`
            If True, searches the Bailer-Jones catalogue and stores its
            parallax equivalent.
        """
        from .utils import search_star, choice_star
        if value not in [True, False]:
            raise AttributeError('bjones attribute must be True or False')
        if value and 'bjones_par' not in self._attributes:
            if hasattr(self, 'code'):
                catalogue = search_star(code=self.code, columns=['**'], catalog='I/347/gaia2dis', verbose=self._verbose)
            else:
                catalogue = search_star(coord=self.coord, columns=['**'], radius=2*u.arcsec,
                                        catalog='I/347/gaia2dis', verbose=self._verbose)
            if len(catalogue) == 0:
                raise ValueError('No star was found in the Bailer-Jones catalogue. It does not exist or VizieR is out.')
            catalogue = catalogue[0]
            if len(catalogue) > 1:
                print('{} stars were found within 2 arcsec from given coordinate.'.format(len(catalogue)))
                print('The list below is sorted by distance. Please select the correct star')
                if hasattr(self.mag, 'G'):
                    print('Star G mag: {}'.format(self.mag['G']))
                catalogue = choice_star(catalogue, self.coord, ['RA_ICRS', 'DE_ICRS', 'Source'], source='bjones')
                if catalogue is None:
                    return
            self._attributes['bjones_par'] = ((1.0/catalogue['rest'][0])*u.arcsec).to(u.mas)
            self.meta_bjones = {c: catalogue[c][0] for c in catalogue.columns}
        self._bjones = value

    @property
    def diameter_gaia(self):
        """Returns the Gaia angular diameter estimate when available."""
        try:
            rad = self.meta_catalogue.get('Rad')
            if rad is not None and not np.ma.core.is_masked(rad):
                return 2*np.arctan((rad*u.solRad)/self.distance).to(u.mas)
        except:
            return None
