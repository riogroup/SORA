from abc import ABCMeta, abstractmethod
from dataclasses import dataclass

import astropy.units as u
import numpy as np
import pyvo
import requests
from astropy.coordinates import SkyCoord
from astropy.time import Time
from astroquery.vizier import Vizier

from sora.config.input_tests import SelectDefault


@dataclass
class Catalogue(metaclass=ABCMeta):
    """Base class describing the columns and services used by a star catalogue."""

    name: str
    cat_path: str
    code: str
    ra: str
    dec: str
    epoch: (str, Time)
    pmra: str = None
    pmdec: str = None
    parallax: str = None
    rad_vel: str = None
    band: dict = None
    errors: list = None
    correlations: dict = None
    ruwe: str = None
    duplicated_source: str = None

    @property
    def _catalogue_service(self):
        """Returns the service used by the catalogue."""
        if hasattr(self, 'tap_url'):
            return 'LIneA TAP'
        return 'VizieR'

    @property
    def _catalogue_reference(self):
        """Returns the catalogue reference in its service."""
        return self.cat_path

    @property
    def _catalogue_description(self):
        """Returns a human-readable catalogue description."""
        return f'{self.name} ({self._catalogue_service}: {self._catalogue_reference})'

    @abstractmethod
    def search_star(self):
        """Searches for one star in the catalogue."""
        pass

    @abstractmethod
    def search_region(self, coord, radius=None, width=None, height=None, columns=None, verbose=False,
                      row_limit=10000000, timeout=600, **kwargs):
        """Searches for stars in a sky region."""
        pass

    def parse_catalogue(self, table):
        """Parses a catalogue table into SORA star parameters.

        Parameters
        ----------
        table : `astropy.table.Table`
            Table with the parameters read from the catalogue server.

        Returns
        -------
        cat_info : `dict`
            Dictionary with catalogue parameters converted to SORA names and units.
        """
        output: dict = {'code': table[self.code].tolist()}
        params = ['ra', 'dec', 'pmra', 'pmdec', 'parallax', 'rad_vel']
        units = [u.deg, u.deg, u.mas / u.year, u.mas / u.year, u.mas, u.km / u.s]
        errors = {}
        if self.errors is None:
            self.errors = [None]*6
        for p, unit, err in zip(params, units, self.errors):
            if not getattr(self, p, None):
                vals = np.zeros(len(table)) * unit
            else:
                vals = self._as_quantity(table, getattr(self, p), unit)
            vals[np.where(np.isnan(vals))] = 0 * unit
            output[p] = vals
            if err is not None and err in table.colnames:
                vals = self._as_quantity(table, err, unit)
                vals[np.where(np.isnan(vals))] = 0 * unit
                errors[p] = vals
        output['errors'] = errors
        if isinstance(self.epoch, Time):
            output['epoch'] = Time(np.repeat(self.epoch, len(table)))
        else:
            output['epoch'] = Time(self._as_array(table, self.epoch), format='jyear')
        if self.band is not None:
            bands = {}
            for keys, item in self.band.items():
                if item not in table.columns:
                    continue
                band_unit = u.Unit(table[item].unit) if table[item].unit is not None else u.mag
                bands[keys] = self._as_quantity(table, item, band_unit)
            output['band'] = bands
        return output

    @staticmethod
    def _column_name(column):
        """Returns a normalized table column name."""
        return column.replace('(', '_').replace(')', '_')

    def _as_array(self, table, column):
        """Returns a table column as a float array with masked values as NaN."""
        values = np.ma.asarray(table[self._column_name(column)])
        if np.ma.isMaskedArray(values):
            values = values.filled(np.nan)
        return np.asarray(values, dtype=float)

    def _as_quantity(self, table, column, unit):
        """Returns a table column as an astropy Quantity."""
        return self._as_array(table, column) * unit

    def get_simple_columns(self):
        """Returns the minimum set of columns required by SORA."""
        col_keys = ['code', 'ra', 'dec', 'pmra', 'pmdec', 'parallax', 'rad_vel', 'epoch']
        columns = [getattr(self, p) for p in col_keys if isinstance(getattr(self, p, None), str)]
        if self.band:
            columns.extend(self.band.values())
        return list(dict.fromkeys(columns))

    def get_choice_columns(self):
        """Returns the columns used when the user needs to choose one source."""
        columns = [self.ra, self.dec]
        if self.band:
            first_band = next(iter(self.band.values()), None)
            if first_band is not None:
                columns.append(first_band)
        return columns


class VizierCatalogue(Catalogue):
    """Defines parameters needed to download star information from VizieR.

    Parameters
    ----------
    name : `str`
        Name of the catalogue used by other processes.
    cat_path : `str`
        The path of the catalogue in the Vizier website.
        For instance, for GaiaEDR3, ``cat_path='I/350/gaiaedr3'``
    code : `str`
        Column name for the unique source identifier within the catalogue.
    ra : `str`
        Column name for right ascension within the catalogue.
    dec : `str`
        Column name for declination within the catalogue.
    epoch : `str`, `astropy.time.Time`
        The epoch of the catalogue. If it is defined in the catalogue, just pass
        the keyword within the catalogue. If the epoch is not present in the catalogue
        table, we must pass a Time object directly, for example ``epoch=Time('J2000')``,
        which defines the catalogue coordinates in J2000 TDB.
    pmra : `str`, optional
        Column name for the proper motion in RA*cos(DEC) within the catalogue.
        If not available, set it to None.
    pmdec : `str`, optional
        Column name for the proper motion in DEC within the catalogue.
        If not available, set it to None.
    parallax : `str`, optional
        Column name for parallax within the catalogue.
        If not available, set it to None.
    rad_vel : `str`, optional
        Column name for radial velocity within the catalogue.
    band : `dict` [`str`, `str`], optional
        A dictionary where the key is band name and the value is the
        keyword referring to the band within the catalogue.
        For instance, in Gaia: ``band={'G': 'Gmag'}``.
        If not available, set it to None.
    errors : `list` [`str`], optional
        A list with the 6 keywords that refer to the uncertainty parameters
        within the catalogue in the order: [ra, dec, pmra, pmdec, parallax, rad_vel].
        If some parameters are not available, please pass each one as None.
        Ex: ``errors=['eRA', 'eDEC', None, None, None, None]``, or ``errors=None`` if
        none of the errors is available.

    Examples
    --------

    To define the Gaia-EDR3 catalogue with `VizierCatalogue`, define an object like:

    >>> catalogue = VizierCatalogue(name='GaiaEDR3', cat_path='I/350/gaiaedr3', code='Source', ra='RA_ICRS', dec='DE_ICRS',
    >>>                             pmra='pmRA', pmdec='pmDE', epoch='Epoch', parallax='Plx', rad_vel='RVDR2', band={'G': 'Gmag'},
    >>>                             errors=['e_RA_ICRS', 'e_DE_ICRS', 'e_pmRA', 'e_pmDE', 'e_Plx', 'e_RVDR2'])

    """

    def __init__(self, **kwargs):
        """Initializes a VizierCatalogue object."""
        super(VizierCatalogue, self).__init__(**kwargs)

    def search_star(self, code=None, coord=None, radius=None):
        """Searches for a specific star in the catalogue.

        Parameters
        ----------
        code : `str`
            Unique source identifier of the star.
        coord : `str`, `astropy.coordinates.SkyCoord`
            Target coordinate to search. It may be specified as a string, in which
            case it is resolved using online services or as the appropriate astropy SkyCoord object.
            ICRS coordinates may also be entered as a string.
        radius : `number`
            Radius of the circular region to query.

        Returns
        -------
        catalogue : `astroquery.utils.commons.TableList`
            Query result with catalogue information about the star.

        Raises
        ------
        ValueError
            Raised when ``code`` or (``coord``, ``radius``) are not provided.

        Notes
        -----
        This function must be called in one of the following ways:

            - Using ``code`` if the unique id of the star is known
            - Using ``coord`` **and** ``radius`` if the catalogue position is known.

            If both alternatives are provided, only the first is used.
        """
        if code is not None:
            vquery = Vizier(columns=['**'], timeout=600)
            kwargs = {self.code: code}
            catalogue = vquery.query_constraints(catalog=self.cat_path, cache=False, **kwargs)
        elif coord is not None:
            catalogue = self.search_region(coord=coord, radius=radius)
        else:
            raise ValueError('At least a code or coord should be given as input')
        # TODO(Implement choice star if necessary)
        return catalogue

    def search_region(self, coord, radius=None, width=None, height=None, columns=None,
                      row_limit=10_000_000, timeout=600, **kwargs):
        """Searches the catalogue around a sky position.

        Parameters
        ----------
        coord : `str`, `astropy.coordinates.SkyCoord`
            Target around which to search. It may be specified as a string, in which
            case it is resolved using online services or as the appropriate astropy SkyCoord object.
            ICRS coordinates may also be entered as a string.
        radius : `number`
            Radius of the circular region to query.
        width : `number`
            Width of the square or rectangular region to query.
        height : `number`
            When set in addition to ``width``, the queried region becomes
            rectangular, with the specified ``width`` and ``height``.
        columns : `list`
            List of strings with the keywords to fetch from the catalogue.
            If ``columns=None`` it will download all the columns.
            If ``columns="simple"`` it will download only the columns for
            the code, epoch and astrometric parameters.
        row_limit : `int`
            Maximum number of rows that will be fetched from the result
            (set to -1 for unlimited). Default: ``row_limit=10_000_000``.
        timeout : `number`
            Timeout for connecting to server in seconds. Default: ``timeout=600``.
        **kwargs
            Additional keyword arguments passed to `astroquery.vizier.Vizier`.

        Returns
        -------
        catalogue : `astroquery.utils.commons.TableList`
            Query result with catalogue information about the stars.
        """
        if columns == 'simple':
            columns = self.get_simple_columns()
        elif columns is None:
            columns = ['**']
        vquery = Vizier(columns=columns, row_limit=row_limit, timeout=timeout, **kwargs)
        catalogue = vquery.query_region(coord, radius=radius, width=width, height=height, catalog=self.cat_path, cache=False)
        return catalogue

    def __repr__(self):
        """Returns the object representation."""
        return self.__str__()

    def __str__(self):
        """Returns the printable catalogue description."""
        return f'<VizierCatalogue: {self.name} defined in https://vizier.cds.unistra.fr/viz-bin/VizieR-3?-source={self.cat_path}>'


class LineaGaiaCatalogue(Catalogue):
    """Gaia catalogue served by the LIneA TAP service."""

    def __init__(self, tap_url='https://userquery.linea.org.br/tap', language=None, **kwargs):
        """Initializes a LineaGaiaCatalogue object.

        Parameters
        ----------
        tap_url : `str`, optional
            URL of the TAP service.
        language : `str`, optional
            Query language passed to the TAP service.
        **kwargs
            Catalogue metadata passed to `Catalogue`.
        """
        self.tap_url = tap_url
        self.language = language
        super(LineaGaiaCatalogue, self).__init__(**kwargs)

    @staticmethod
    def _format_value(value):
        """Formats a value for use in a TAP query."""
        text = str(value).strip()
        if text.isdigit():
            return text
        return "'{}'".format(text.replace("'", "''"))

    def _run_query(self, query):
        """Runs a TAP query and returns the resulting table."""
        session = requests.Session()
        tap = pyvo.dal.TAPService(self.tap_url, session=session)
        if self.language is None:
            return tap.run_sync(query).to_table()
        return tap.run_sync(query, language=self.language).to_table()

    @staticmethod
    def _format_float(value):
        """Formats a floating-point value for use in ADQL."""
        return f'{float(value):.16f}'

    def _format_columns(self, columns):
        """Formats selected columns for a TAP SELECT clause."""
        if columns == 'simple':
            columns = self.get_simple_columns()
        if columns is None:
            return '*'
        return ', '.join(dict.fromkeys(columns))

    def _format_column_filters(self, column_filters):
        """Formats column filters for a TAP WHERE clause."""
        if not column_filters:
            return []
        filters = []
        for column, expression in column_filters.items():
            filters.append(f'{column} {str(expression).strip()}')
        return filters

    def search_star(self, code=None, coord=None, radius=None):
        """Searches for a specific star in the LIneA Gaia catalogue.

        Parameters
        ----------
        code : `str`
            Unique source identifier of the star.
        coord : `str`, `astropy.coordinates.SkyCoord`
            Target coordinate used when searching by sky position.
        radius : `astropy.units.Quantity`, optional
            Radius of the circular region to query when ``coord`` is used.

        Returns
        -------
        catalogue : `list`
            Empty list when no source is found, or a list containing one table.
        """
        if code is not None:
            query = f'SELECT * FROM {self.cat_path} WHERE {self.code} = {self._format_value(code)}'
            catalogue = self._run_query(query)
        elif coord is not None:
            catalogue = self.search_region(coord=coord, radius=radius)
        else:
            raise ValueError('At least a code or coord should be given as input')
        if isinstance(catalogue, list):
            return catalogue
        return [] if len(catalogue) == 0 else [catalogue]

    def search_region(self, coord, radius=None, width=None, height=None, columns=None,
                      row_limit=10_000_000, timeout=600, column_filters=None, **kwargs):
        """Searches the LIneA Gaia catalogue around a sky position.

        Parameters
        ----------
        coord : `str`, `astropy.coordinates.SkyCoord`
            Target around which to search.
        radius : `astropy.units.Quantity`, optional
            Radius of the circular region to query.
        width : `astropy.units.Quantity`, optional
            Width of the rectangular region to query.
        height : `astropy.units.Quantity`, optional
            Height of the rectangular region to query. If omitted, ``width`` is
            used for both dimensions.
        columns : `list`, `str`, optional
            Columns to retrieve. Use ``'simple'`` for the minimum SORA columns.
        row_limit : `int`, optional
            Maximum number of rows to fetch. Non-positive values remove the
            ``TOP`` clause.
        timeout : `number`, optional
            Accepted for API compatibility with VizieR catalogues.
        column_filters : `dict`, optional
            Additional ADQL filter expressions keyed by column name.
        **kwargs
            Accepted for API compatibility with VizieR catalogues.

        Returns
        -------
        catalogue : `list`
            Empty list when no source is found, or a list containing one table.
        """
        del timeout, kwargs

        if not isinstance(coord, SkyCoord):
            coord = SkyCoord(coord, unit=('hourangle', 'deg'))
        coord = coord.icrs
        ra_deg = self._format_float(coord.ra.deg)
        dec_deg = self._format_float(coord.dec.deg)

        if radius is not None:
            radius_deg = self._format_float(u.Quantity(radius).to_value(u.deg))
            region_filter = (
                f"1=CONTAINS(POINT('ICRS', {self.ra}, {self.dec}), "
                f"CIRCLE('ICRS', {ra_deg}, {dec_deg}, {radius_deg}))"
            )
        elif width is not None:
            width_deg = self._format_float(u.Quantity(width).to_value(u.deg))
            height_deg = self._format_float(u.Quantity(height if height is not None else width).to_value(u.deg))
            region_filter = (
                f"1=CONTAINS(POINT('ICRS', {self.ra}, {self.dec}), "
                f"BOX('ICRS', {ra_deg}, {dec_deg}, {width_deg}, {height_deg}))"
            )
        else:
            raise ValueError('At least a radius or width should be given as input')

        query_filters = [region_filter]
        query_filters.extend(self._format_column_filters(column_filters))
        top = f'TOP {int(row_limit)} ' if row_limit and row_limit > 0 else ''
        query = f"SELECT {top}{self._format_columns(columns)} FROM {self.cat_path} WHERE {' AND '.join(query_filters)}"
        catalogue = self._run_query(query)
        return [] if len(catalogue) == 0 else [catalogue]

    def __repr__(self):
        """Returns the object representation."""
        return self.__str__()

    def __str__(self):
        """Returns the printable catalogue description."""
        return f'<LineaGaiaCatalogue: {self.name} defined in {self.tap_url} ({self.cat_path})>'


def is_timeout_error(exc):
    """Returns True when an exception indicates a timeout condition."""
    if isinstance(exc, (requests.exceptions.Timeout, TimeoutError)):
        return True
    text = str(exc).lower()
    return 'timed out' in text or 'timeout' in text


def should_fallback_to_gaiadr3(catalogue, exc):
    """Returns True when a LIneA TAP query should fall back to VizieR Gaia DR3."""
    return isinstance(catalogue, LineaGaiaCatalogue) and is_timeout_error(exc)


gaiadr2 = VizierCatalogue(name='GaiaDR2', cat_path='I/345/gaia2', code='Source', ra='RA_ICRS', dec='DE_ICRS',
                          pmra='pmRA', pmdec='pmDE', epoch='Epoch', parallax='Plx', rad_vel='RV', band={'G': 'Gmag'},
                          errors=['e_RA_ICRS', 'e_DE_ICRS', 'e_pmRA', 'e_pmDE', 'e_Plx', 'e_RV'],
                          correlations={('ra', 'dec'): 'RADEcor', ('ra', 'parallax'): 'RAPlxcor',
                                        ('ra', 'pmra'): 'RApmRAcor', ('ra', 'pmdec'): 'RApmDEcor',
                                        ('dec', 'parallax'): 'DEPlxcor', ('dec', 'pmra'): 'DEpmRAcor',
                                        ('dec', 'pmdec'): 'DEpmDEcor', ('parallax', 'pmra'): 'PlxpmRAcor',
                                        ('parallax', 'pmdec'): 'PlxpmDEcor', ('pmra', 'pmdec'): 'pmRApmDEcor'},
                          ruwe='RUWE', duplicated_source='Dup')

gaiaedr3 = VizierCatalogue(name='GaiaEDR3', cat_path='I/350/gaiaedr3', code='Source', ra='RA_ICRS', dec='DE_ICRS',
                           pmra='pmRA', pmdec='pmDE', epoch='Epoch', parallax='Plx', rad_vel='RVDR2', band={'G': 'Gmag'},
                           errors=['e_RA_ICRS', 'e_DE_ICRS', 'e_pmRA', 'e_pmDE', 'e_Plx', 'e_RVDR2'],
                           correlations={('ra', 'dec'): 'RADEcor', ('ra', 'parallax'): 'RAPlxcor',
                                         ('ra', 'pmra'): 'RApmRAcor', ('ra', 'pmdec'): 'RApmDEcor',
                                         ('dec', 'parallax'): 'DEPlxcor', ('dec', 'pmra'): 'DEpmRAcor',
                                         ('dec', 'pmdec'): 'DEpmDEcor', ('parallax', 'pmra'): 'PlxpmRAcor',
                                         ('parallax', 'pmdec'): 'PlxpmDEcor', ('pmra', 'pmdec'): 'pmRApmDEcor'},
                           ruwe='RUWE', duplicated_source='Dup')

epoch = Time('J2016.0', scale='tdb')
gaiadr3 = VizierCatalogue(name='GaiaDR3', cat_path='I/355/gaiadr3', code='Source', ra='RA_ICRS', dec='DE_ICRS',
                          pmra='pmRA', pmdec='pmDE', epoch=epoch, parallax='Plx', rad_vel='RV', band={'G': 'Gmag'},
                          errors=['e_RA_ICRS', 'e_DE_ICRS', 'e_pmRA', 'e_pmDE', 'e_Plx', 'e_RV'],
                          correlations={('ra', 'dec'): 'RADEcor', ('ra', 'parallax'): 'RAPlxcor',
                                        ('ra', 'pmra'): 'RApmRAcor', ('ra', 'pmdec'): 'RApmDEcor',
                                        ('dec', 'parallax'): 'DEPlxcor', ('dec', 'pmra'): 'DEpmRAcor',
                                        ('dec', 'pmdec'): 'DEpmDEcor', ('parallax', 'pmra'): 'PlxpmRAcor',
                                        ('parallax', 'pmdec'): 'PlxpmDEcor', ('pmra', 'pmdec'): 'pmRApmDEcor'},
                          ruwe='RUWE', duplicated_source='Dup')

gaiadr3_linea = LineaGaiaCatalogue(name='GaiaDR3', cat_path='gaia_dr3.source', code='source_id', ra='ra', dec='dec',
                              pmra='pmra', pmdec='pmdec', epoch='ref_epoch', parallax='parallax',
                              rad_vel='radial_velocity',
                              band={'G': 'phot_g_mean_mag'},
                              errors=['ra_error', 'dec_error', 'pmra_error', 'pmdec_error',
                                      'parallax_error', 'radial_velocity_error'],
                              correlations={('ra', 'dec'): 'ra_dec_corr', ('ra', 'parallax'): 'ra_parallax_corr',
                                            ('ra', 'pmra'): 'ra_pmra_corr', ('ra', 'pmdec'): 'ra_pmdec_corr',
                                            ('dec', 'parallax'): 'dec_parallax_corr',
                                            ('dec', 'pmra'): 'dec_pmra_corr', ('dec', 'pmdec'): 'dec_pmdec_corr',
                                            ('parallax', 'pmra'): 'parallax_pmra_corr',
                                            ('parallax', 'pmdec'): 'parallax_pmdec_corr',
                                            ('pmra', 'pmdec'): 'pmra_pmdec_corr'},
                              ruwe='ruwe', duplicated_source='duplicated_source')

epoch_ubsc = Time('J1991.25', scale='tdb')
ubsc = VizierCatalogue(name='UBSC', cat_path='J/AJ/164/36/table2', code='HIP', ra='RA_ICRS', dec='DE_ICRS',
                       pmra='pmRA', pmdec='pmDE', epoch=epoch_ubsc, parallax='plx', band={'Hp': 'Hpmag'},
                       errors=['e_RA_ICRS', 'e_DE_ICRS', 'e_pmRA', 'e_pmDE', 'e_plx', None])

allowed_catalogues = SelectDefault(instance=Catalogue,
                                   defaults={'gaiadr3_linea': gaiadr3_linea, 
                                             'gaiadr2': gaiadr2, 
                                             'gaiaedr3': gaiaedr3,
                                             'gaiadr3': gaiadr3, 
                                             'ubsc': ubsc})
