import warnings

import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.time import Time

from sora.config import input_tests
from sora.config.decorators import deprecated_alias, deprecated_function
from .meta import MetaStar
from .utils import search_star, van_belle, kervella, spatial_motion, choice_star

warnings.simplefilter('always', UserWarning)

__all__ = ['Star']


class Star(MetaStar):
    """Defines a star.

    Parameters
    ----------
    catalogue : `str`
        The catalogue to download data. It can be ``'gaiadr2'`` or ``'gaiaedr3'``.

    code : `str`
        Gaia Source code for searching in VizieR.

    coord : `str`, `astropy.coordinates.SkyCoord`
        If code is not given, coord nust have the coordinates RA and DEC of the
        star to search in VizieR: ``'hh mm ss.ss +dd mm ss.ss'``.

    ra : `int`, `float`
        Right Ascension, in deg.

    dec : `int`, `float`
        Declination, in deg.

    parallax : `int`, `float`. default=0
        Parallax, in mas.

    pmra : `int`, `float`, default=0
        Proper Motion in RA*, in mas/year.

    pmdec : `int`, `float`, default=0
        Proper Motion in DEC, in mas/year.

    rad_vel : `int`, `float`, default=0
        Radial Velocity, in km/s.

    epoch : `str`, `astropy.time.Time`, default='J2000'
        Epoch of the coordinates.

    nomad : `bool`
        If True, it tries to download the magnitudes from NOMAD catalogue.

    bjones : `bool`, default=True
        If True, it uses de star distance from Bailer-Jones et al. (2018).

    verbose : `bool`, default=True
        If True, it prints the downloaded information

    local : `bool`, default=False
        If True, it uses the given coordinate in 'coord' as final coordinate.

    Note
    ----
    The user can give either 'coord' or 'ra' and 'dec', but not both.

    To download the coordinates from Gaia, "local" must be set as False
    and the ("code") or ("coord") or ("ra" and "dec") must be given.

    All values downloaded from Gaia will replace the ones given by the user.

    """

    @deprecated_alias(log='verbose')  # remove this line in v1.0
    def __init__(self, catalogue='gaiaedr3', **kwargs):

        self._attributes = {}
        self.mag = {}
        self.errors = {'RA': 0*u.mas, 'DEC': 0*u.mas, 'Plx': 0*u.mas, 'pmRA': 0*u.mas/u.year,
                       'pmDE': 0*u.mas/u.year, 'rad_vel': 0*u.km/u.year}
        allowed_kwargs = ['bjones', 'code', 'coord', 'dec', 'epoch', 'local', 'verbose', 'nomad', 'parallax', 'pmdec', 'pmra',
                          'ra', 'rad_vel']
        input_tests.check_kwargs(kwargs, allowed_kwargs=allowed_kwargs)
        allowed_catalogues = ['gaiadr2', 'gaiaedr3']
        if catalogue not in allowed_catalogues:
            raise ValueError('Catalogue {} is not one of the allowed catalogues {}'.format(catalogue, allowed_catalogues))
        self._catalogue = {'gaiadr2': 'Gaia-DR2', 'gaiaedr3': 'Gaia-EDR3'}[catalogue]
        self._verbose = kwargs.get('verbose', True)
        local = kwargs.get('local', False)
        self.bjones = False
        if 'code' in kwargs:
            self.code = kwargs['code']
        if 'coord' in kwargs:
            if 'ra' in kwargs or 'dec' in kwargs:
                raise ValueError("User must give 'coord' or 'ra' and 'dec', not both")
            coord = SkyCoord(kwargs['coord'], unit=('hourangle', 'deg'))
            self.ra = coord.ra
            self.dec = coord.dec
        if 'ra' in kwargs and 'dec' in kwargs:
            self.ra = kwargs.get('ra')
            self.dec = kwargs.get('dec')
        self.parallax = kwargs.get('parallax', 0.0)
        self.pmra = kwargs.get('pmra', 0.0)
        self.pmdec = kwargs.get('pmdec', 0.0)
        self.rad_vel = kwargs.get('rad_vel', 0.0)
        self.epoch = kwargs.get('epoch', 'J2000')
        if local:
            if 'RA' not in self._attributes or 'DEC' not in self._attributes:
                raise ValueError("User must give 'ra' and 'dec' for local coordinates")
        else:
            if not hasattr(self, 'code') and 'RA' not in self._attributes:
                raise ValueError("User must give gaia Source ID 'code' or coordinates for the online search")
            self.__searchgaia(catalog=catalogue)
        if kwargs.get('nomad', True):
            self.__getcolors()
        try:
            self.bjones = kwargs.get('bjones', False)
        except ValueError:
            pass

    def set_magnitude(self, **kwargs):
        """Sets the magnitudes of a star.

        Parameters
        ----------
        band=value : `str`
            The star magnitude for given band. The band name can be any string
            the user wants.

        Examples
        --------
        To set the stars magnitude in the band G:\n
        >>> set_magnitude(G=10)

        To set the star's magnitude in the band K:\n
        >>> set_magnitude(K=15)

        To set the star's magnitude in a customized band:\n
        >>> set_magnitude(newband=6)
        """
        for key in kwargs:
            mag = input_tests.test_attr(kwargs[key], float, key)
            if key in self.mag:
                warnings.warn('{0} mag already defined. {0}={1} will be replaced by {0}={2}'.format(
                    key, self.mag[key], mag))
            self.mag[key] = mag

    def set_diameter(self, diameter):
        """Sets an user diameter for the star, in mas.

        Parameters
        ----------
        diameter : `int`, `float`
            Sets the user diameter of the star, in mas.
        """
        self.diameter_user = diameter * u.mas
        if diameter < 0:
            warnings.warn("negative sizes are converted to positive.")
            self.diameter_user = np.absolute(self.diameter_user)

    def van_belle(self):
        """Determines the diameter of a star in mas using equations from van Belle (1999).

        See: Publi. Astron. Soc. Pacific 111, 1515-1523:.
        """
        return van_belle(self.mag.get('B'), self.mag.get('V'), self.mag.get('K'))

    def kervella(self):
        """Determines the diameter of a star in mas using equations from Kervella et. al (2004).

        See: A&A Vol.  426, No.  1:.
        """
        return kervella(self.mag.get('B'), self.mag.get('V'), self.mag.get('K'))

    @deprecated_alias(log='verbose')  # remove this line in v1.0
    def apparent_diameter(self, distance, mode='auto', band='V', star_type='sg', verbose=True):
        """Calculates the apparent diameter of the star at a given distance.

        Parameters
        ----------
        distance : `int`, `float`
            Object geocentric distance, in AU.

        mode : `str`, default='auto'
            The mode to calculate the apparent diameter.\n
            - ``'user'``: calculates using user given diameter.\n
            - ``'gaia'``: calculates using diameter obtained from Gaia.\n
            - ``'kervella'``: calculates using Kervella equations.\n
            - ``'van_belle'``: calculates using van Belle equations.\n
            - ``'auto'``: tries all the above methods until it is able to calculate diameter.\n
            The order of try is the same as shown above (user, Gaia, Kervella, Van Belle).

        band : `str`
            The band filter to calculate the diameter. If mode is `kervella`
            or `van_belle`, the filter must be given, ``'B'`` or ``'V'``.
            If mode `auto`, ``'V'`` is selected.

        star_type :`str`
            Type of star to calculate the diameter. If mode is `van_belle`,
            the star type must be given. If mode is `auto`, ``star_type='sg'``.\n
            Accepted types:\n
            - ``'sg'`` for 'Super Giant'.\n
            - ``'ms'`` for 'Main Sequence'.\n
            - ``'vs'`` for 'Variable Star'.

        verbose : `bool`
            If True, it prints the mode used by `auto`.
        """
        try:
            distance = distance.to(u.km)
        except:
            distance = distance * u.AU
        if distance < 0:
            warnings.warn("negative distances are converted to positive.")
            distance = np.absolute(distance)

        if mode in ['user', 'auto']:
            try:
                diam = distance*np.tan(self.diameter_user)
                if verbose:
                    print('Calculating apparent diameter from user defined diameter')
                return diam.to(u.km)
            except:
                pass

        if mode == 'user':
            raise ValueError('User diameter must be informed.')

        if mode in ['gaia', 'auto']:
            try:
                diam = distance*np.tan(self.diameter_gaia)
                if verbose:
                    text = ''
                    if self.bjones:
                        text += ' + Bailer-Jones et al. (2018)'
                    print('Apparent diameter using Gaia' + text)
                return diam.to(u.km)
            except:
                pass

        if mode == 'gaia':
            raise ValueError('It is not possible to calculate star diameter from Gaia.')

        if band not in ['B', 'V']:
            raise KeyError('band must be informed as "B", or "V"')

        if mode in ['kervella', 'auto']:
            diam_kerv = self.kervella().get(band)
            if diam_kerv is None:
                raise ValueError('Diameter could not be calculated for given band')
            if verbose:
                print('Apparent diameter using Kervella et al. (2004)')
            diam = distance*np.tan(diam_kerv)
            return diam.to(u.km)

        if star_type not in ['sg', 'ms', 'vs']:
            raise KeyError('star_type must be informed as "sg", "ms" or "vs"')

        if mode in ['van_belle', 'auto']:
            diam_van = self.van_belle().get(star_type)
            if diam_van is None:
                raise ValueError('Diameter could not be calculated using Van Belle')
            diam_van = diam_van.get(band)
            if diam_van is None:
                raise ValueError('Diameter could not be calculated for given band')
            if verbose:
                print('Apparent diameter using van Belle (1999)')
            diam = distance*np.tan(diam_van)
            return diam.to(u.km)

        raise AttributeError("Star apparent diameter could not be calculated. ",
                             "Please define star diameter or B,V,K magnitudes.")

    def __searchgaia(self, catalog):
        """Searches for the star position in the Gaia catalogue and save information.

        Parameters
        ----------
        catalog : `str`
            The catalogue to download data. It can be ``'gaiadr2'`` or ``'gaiaedr3'``.
        """
        catalogues = {'gaiadr2': 'I/345/gaia2', 'gaiaedr3': 'I/350/gaiaedr3'}
        cat = catalogues[catalog]
        if hasattr(self, 'code'):
            catalogue = search_star(code=self.code, columns=['**'], catalog=cat, verbose=self._verbose)
        else:
            catalogue = search_star(coord=self.coord, columns=['**'], radius=2*u.arcsec,
                                    catalog=cat, verbose=self._verbose)
        if len(catalogue) == 0:
            raise ValueError('No star was found. It does not exist or VizieR is out.')
        catalogue = catalogue[0]
        if len(catalogue) > 1:
            if self._verbose:
                print('{} stars were found within 2 arcsec from given coordinate.'.format(len(catalogue)))
                print('The list below is sorted by distance. Please select the correct star')
            catalogue = choice_star(catalogue, self.coord, ['RA_ICRS', 'DE_ICRS', 'Gmag'], source='gaia')
        self.code = catalogue['Source'][0]
        self.ra = catalogue['RA_ICRS'][0]*u.deg
        self.dec = catalogue['DE_ICRS'][0]*u.deg
        self.pmra = catalogue['pmRA'][0]*u.mas/u.year
        self.pmdec = catalogue['pmDE'][0]*u.mas/u.year
        self.epoch = Time(catalogue['Epoch'][0], format='jyear')
        self.parallax = catalogue['Plx'][0]*u.mas
        rv_name = {'gaiadr2': 'RV', 'gaiaedr3': 'RVDR2'}
        self.rad_vel = catalogue[rv_name[catalog]][0]*u.km/u.s
        self.set_magnitude(G=catalogue['Gmag'][0])

        self.meta_gaia = {c: catalogue[c][0] for c in catalogue.columns}

        self.errors['RA'] = self.meta_gaia['e_RA_ICRS']*u.mas
        self.errors['DEC'] = self.meta_gaia['e_DE_ICRS']*u.mas
        self.errors['Plx'] = self.meta_gaia['e_Plx']*u.mas
        self.errors['pmRA'] = self.meta_gaia['e_pmRA']*(u.mas/u.yr)
        self.errors['pmDE'] = self.meta_gaia['e_pmDE']*(u.mas/u.yr)
        erv_name = {'gaiadr2': 'e_RV', 'gaiaedr3': 'e_RVDR2'}
        self.errors['rad_vel'] = self.meta_gaia[erv_name[catalog]]*(u.km/u.s)

        A = (1*u.AU).to(u.km).value
        cov = np.zeros((6, 6))
        a = ['RA', 'DE', 'Plx', 'pmRA', 'pmDE']
        for i in np.arange(5):
            v1 = 'e_' + a[i]
            if i in [0, 1]:
                v1 += '_ICRS'
            for j in np.arange(i, 5):
                v2 = 'e_' + a[j]
                if j in [0, 1]:
                    v2 += '_ICRS'
                if i == j:
                    x = self.meta_gaia[v1]**2
                    if not np.ma.core.is_masked(x):
                        cov[i, i] = x
                else:
                    x = self.meta_gaia[a[i]+a[j]+'cor']*self.meta_gaia[v1]*self.meta_gaia[v2]
                    if not np.ma.core.is_masked(x):
                        cov[i, j] = x
                        cov[j, i] = cov[i, j]
            x = cov[i, 2]*(self.meta_gaia[rv_name[catalog]]/A)
            if not np.ma.core.is_masked(x):
                cov[i, 5] = x
                cov[5, i] = cov[i, 5]
        x = cov[2, 2]*(self.meta_gaia[rv_name[catalog]]**2 + self.meta_gaia[erv_name[catalog]]**2)/(A**2) \
            + (self.meta_gaia['Plx']*self.meta_gaia[erv_name[catalog]]/A)**2
        if not np.ma.core.is_masked(x):
            cov[5, 5] = x
        cov[np.where(np.isnan(cov))] = 0.0
        self.cov = cov

        if self._verbose:
            print('1 {} star found G={}'.format(self._catalogue, catalogue['Gmag'][0]))
            print('star coordinate at J{}: RA={} +/- {}, DEC={} +/- {}'.format(self.epoch.jyear,
                  self.ra.to_string(u.hourangle, sep='hms', precision=5), self.errors['RA'],
                  self.dec.to_string(u.deg, sep='dms', precision=4), self.errors['DEC']))

    def __getcolors(self):
        """ Searches for the B,V,K magnitudes of the star in the NOMAD catalogue on VizieR.
        """
        columns = ['RAJ2000', 'DEJ2000', 'Bmag', 'Vmag', 'Rmag', 'Jmag', 'Hmag', 'Kmag']
        catalogue = search_star(coord=self.coord, columns=columns, radius=2*u.arcsec,
                                catalog='I/297/out', verbose=self._verbose)
        if len(catalogue) == 0:
            if self._verbose:
                warnings.warn('No star was found on NOMAD that matches the star')
            return
        catalogue = catalogue[0]
        if len(catalogue) > 1:
            print('{} stars were found within 2 arcsec from given coordinate.'.format(len(catalogue)))
            print('The list below is sorted by distance. Please select the correct star')
            if hasattr(self.mag, 'G'):
                print('Star G mag: {}'.format(self.mag['G']))
            catalogue = choice_star(catalogue, self.coord, ['RAJ2000', 'DEJ2000', 'Bmag', 'Vmag',
                                                            'Rmag', 'Jmag', 'Hmag', 'Kmag'], source='nomad')
            if catalogue is None:
                return
        errors = []
        for mag in ['B', 'V', 'R', 'J', 'H', 'K']:
            name = mag + 'mag'
            if np.ma.core.is_masked(catalogue[name][0]):
                errors.append(mag)
                continue
            self.set_magnitude(**{mag: catalogue[name][0]})
        if len(errors) > 0 and self._verbose:
            print('Magnitudes in {} were not located in NOMAD'.format(errors))

    # remove this block for v1.0
    @deprecated_function(message="Please use get_position(time=time, observer='geocenter')")
    def geocentric(self, time):
        """ Calculates the position of the star, propagating the position using parallax and proper motion

        Parameters
        ----------
        time : `str`, `astropy.time.Time`
            Reference time to apply proper motion and calculate parallax. It can be a string
            in the ISO format (yyyy-mm-dd hh:mm:ss.s) or an astropy Time object.
        """
        return self.get_position(time=time, observer='geocenter')

    @deprecated_function(message="Please use get_position(time=time, observer='barycenter')")
    def barycentric(self, time):
        """Calculates the position of the star using proper motion.

        Parameters
        ----------
        time : `str`, `astropy.time.Time`
            Reference time to apply proper motion. It can be a string
            in the ISO format (yyyy-mm-dd hh:mm:ss.s) or an astropy Time object.
        """
        return self.get_position(time=time, observer='barycenter')
    # end of block removal

    def get_position(self, time, observer='geocenter'):
        """Calculates the position of the star for given observer,
        propagating the position using parallax and proper motion

        Parameters
        ----------
        time : `float`, `astropy.time.Time`
            Reference time to apply proper motion and calculate parallax. It can be a string
            in the ISO format (yyyy-mm-dd hh:mm:ss.s) or an astropy Time object.

        observer : `str`, `sora.observer.Observer`, `sora.observer.Spacecraft`
            Observer of the star t calculate position. It can be 'geocenter' for a geocentric
            coordinate, 'barycenter' for a barycenter coordinate, or a sora observer object.

        Returns
        -------
        coord : `astropy.coordinates.SkyCoord`
            Astropy SkyCoord object with the star coordinates at the given time.
        """
        from astropy.coordinates import SphericalRepresentation, SkyOffsetFrame, ICRS
        from sora import Observer, Spacecraft

        try:
            time = Time(time)
        except:
            time = Time(time, format='jd', scale='utc')
        if observer not in ['geocenter', 'barycenter'] and not isinstance(observer, (Observer, Spacecraft)):
            raise ValueError("'observer' must be an Observer object or one of the following"
                             " strings: ['geocenter', 'barycenter]")

        def apply_offset(coord):
            if not hasattr(self, 'offset'):
                return coord
            star_frame = SkyOffsetFrame(origin=coord)
            new_pos = SkyCoord(lon=self.offset.d_lon_coslat, lat=self.offset.d_lat, frame=star_frame)
            p = new_pos.transform_to(ICRS)
            return SkyCoord(ra=p.ra, dec=p.dec, distance=p.distance)

        dt = time - self.epoch
        if not time.isscalar:
            if time.max() - time.min() > 1*u.day:
                raise ValueError('list of times must be in a interval of 1 day to process.')
            dt = dt[0]
        bar_star = spatial_motion(self.ra, self.dec, self.pmra, self.pmdec, self.parallax, self.rad_vel, dt=dt.jd)

        if observer == "barycenter" or self.coord.distance.unit.is_unity() or np.isnan(self.coord.distance):
            return apply_offset(bar_star)

        if observer == "geocenter":
            observer = Observer(code='500', ephem='horizons')

        bar_obs = observer.get_vector(time=time, origin='barycenter')

        topo = bar_star.cartesian - bar_obs.cartesian
        topo = topo.represent_as(SphericalRepresentation)
        topo = SkyCoord(topo.lon, topo.lat, topo.distance)
        return apply_offset(topo)

    def error_at(self, time):
        """Estimates the star position error at a given time.

        Parameters
        ----------
        time : `str`, `astropy.time.Time`
            Reference time to project star error. It can be a string
            in the ISO format (yyyy-mm-dd hh:mm:ss.s) or an astropy Time object.

        Returns
        -------
        errors : `list`
            In RA* and DEC.

        """
        try:
            time = Time(time)
        except:
            time = Time(time, format='jd', scale='utc')
        dt = time - self.epoch
        n_coord, errors = spatial_motion(self.ra, self.dec, self.pmra, self.pmdec, self.parallax,
                                         self.rad_vel, dt=dt.jd, cov_matrix=self.cov)
        return errors[0]*u.mas, errors[1]*u.mas

    def add_offset(self, da_cosdec, ddec):
        """Adds an offset to the star position.

        Parameters
        ----------
        da_cosdec : `int`, `float`
            Offset in Delta_alpha_cos_delta, in mas.

        ddec : `int`, `float`
            Offset in Delta_delta, in mas.
        """
        from astropy.coordinates import SphericalCosLatDifferential

        dadc = input_tests.test_attr(da_cosdec, float, 'da_cosdec')
        dd = input_tests.test_attr(ddec, float, 'ddec')
        self.offset = SphericalCosLatDifferential(dadc * u.mas, dd * u.mas, 0.0 * u.km)

    def __str__(self):
        """String representation of the Star class.
        """
        out = ''
        if hasattr(self, 'code'):
            out += '{} star Source ID: {}\n'.format(self._catalogue, self.code)
        else:
            out += 'User coordinates\n'
        out += ('ICRS star coordinate at J{}:\n'
                'RA={} +/- {:.4f}, DEC={} +/- {:.4f}\n'
                'pmRA={:.3f} +/- {:.3f} mas/yr, pmDEC={:.3f} +/- {:.3f} mas/yr\n'
                'Plx={:.4f} +/- {:.4f} mas, Rad. Vel.={:.2f} +/- {:.2f} km/s \n\n'.format(
                    self.epoch.jyear, self.ra.to_string(u.hourangle, sep='hms', precision=5),
                    self.errors['RA'], self.dec.to_string(u.deg, sep='dms', precision=4), self.errors['DEC'],
                    self.pmra.value, self.errors['pmRA'].value, self.pmdec.value, self.errors['pmDE'].value,
                    self.parallax.value, self.errors['Plx'].value, self.rad_vel.value, self.errors['rad_vel'].value))
        if hasattr(self, 'offset'):
            out += 'Offset Apllied: d_alpha_cos_dec = {}, d_dec = {}\n'.format(
                self.offset.d_lon_coslat, self.offset.d_lat)
        out += 'Magnitudes:'
        mag_out = [' {}: {:6.3f}'.format(mag, self.mag[mag]) for mag in self.mag]
        out_mag = []
        for i, mag in enumerate(mag_out):
            if i % 6 == 0:
                out_mag.append([])
            out_mag[-1].append(mag)
        out += (',\n'+' '*11).join([','.join(out_i) for out_i in out_mag])
        out += '\n\n'
        if self.diameter_gaia is not None:
            text = ''
            if self.bjones:
                text += ' + Bailer-Jones et al. (2018)'
            out += 'Apparent diameter: {:.4f}, Source: Gaia-DR2{}\n'.format(self.diameter_gaia, text)
        if hasattr(self, 'diameter_user'):
            out += 'Apparent diameter: {:.4f}, Source: User\n'.format(self.diameter_user)
        kerv = self.kervella()
        if kerv:
            out += 'Apparent diameter from Kervella et. al (2004):\n'
            out += '   ' + ','.join([' {}: {:.4f}'.format(k, v) for k, v in kerv.items()])
        vanb = self.van_belle()
        if vanb:
            out += '\nApparent diameter from van Belle (1999):'
            for key, value in vanb.items():
                out += '\n    {}:'.format(key)
                out += ','.join([' {}: {:.4f}'.format(k, v) for k, v in value.items()])
        return out
