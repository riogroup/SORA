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

    @abstractmethod
    def search_star(self):
        pass

    @abstractmethod
    def search_region(self, coord, radius=None, width=None, height=None, columns=None, verbose=False,
                      row_limit=10000000, timeout=600, **kwargs):
        pass

    def parse_catalogue(self, table):
        """Properly interprets the table

        Parameters
        ----------
        table : `astropy.table.Table`
            The table with the parameters read from the catalogue server

        Returns
        -------
        cat_info : `dict`
            Dictionary with the list of parameters
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
        return column.replace('(', '_').replace(')', '_')

    def _as_array(self, table, column):
        values = np.ma.asarray(table[self._column_name(column)])
        if np.ma.isMaskedArray(values):
            values = values.filled(np.nan)
        return np.asarray(values, dtype=float)

    def _as_quantity(self, table, column, unit):
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
    """VizierCatalogue defines the parameters necessary to download all
    the information of stars from a catalogue on the Vizier webservice.

    Parameters
    ----------
    name : `str`
        The name of the catalogue which will be referred in other processes
    cat_path : `str`
        The path of the catalogue in the Vizier website.
        For instance, for GaiaEDR3, ``cat_path='I/350/gaiaedr3'``
    code : `str`
        The keyword referring to a unique code within the catalogue
    ra : `str`
        The keyword referring to the Right Ascension within the catalogue
    dec : `str`
        The keyword referring to the Declination within the catalogue
    epoch : `str`, `astropy.time.Time`
        The epoch of the catalogue. If it is defined in the catalogue, just pass
        the keyword within the catalogue. If the epoch is not present in the catalogue
        table, we must pass a Time object directly, for example ``epoch=Time('J2000')``,
        which defines the catalogue coordinates in J2000 TDB.
    pmra : `str`, optional
        The keyword referring to the Proper Motion in ra*cosdec within the catalogue.
        If not available, set it to None.
    pmdec : `str`, optional
        The keyword referring to the Proper Motion in dec within the catalogue.
        If not available, set it to None.
    parallax : `str`, optional
        The keyword referring to the Parallax within the catalogue.
        If not available, set it to None.
    rad_vel : `str`, optional
        The keyword referring to the Radial Velocity within the catalogue
    bands : `dict` [`str`, `str`], optional
        A dictionary where the key is band name and the value is the
        keyword referring to the band within the catalogue.
        For instance, in Gaia: ``bands={'G': 'Gmag'}``.
        If not available, set it to None.
    errors : `list` [`str`], optional
        A list with the 6 keywords that refer to the uncertainty parameters
        within the catalogue in the order: [ra, dec, pmra, pmdec, parallax, rad_vel].
        If some parameters are not available, please pass each one as None.
        Ex: ``errors=['eRA', 'eDEC', None, None, None, None]``, or ``errors=None`` if
        none of the errors is available.

    Examples
    --------

    To define the Gaia-EDR3 catalogue with VizierCatalogue, we must define a object like:

    >>> catalogue = VizierCatalogue(name='GaiaEDR3', cat_path='I/350/gaiaedr3', code='Source', ra='RA_ICRS', dec='DE_ICRS',
    >>>                             pmra='pmRA', pmdec='pmDE', epoch='Epoch', parallax='Plx', rad_vel='RVDR2', band={'G': 'Gmag'},
    >>>                             errors=['e_RA_ICRS', 'e_DE_ICRS', 'e_pmRA', 'e_pmDE', 'e_Plx', 'e_RVDR2'])

    """

    def __init__(self, **kwargs):
        """"""
        super(VizierCatalogue, self).__init__(**kwargs)

    def search_star(self, code=None, coord=None, radius=None):
        """Looks for a specific star in the catalogue

        Parameters
        ----------
        code : `str`
            The unique id of the star
        coord : `str`, `astropy.coordinates.SkyCoord`
            The target coordinate which to search. It may be specified as a string in which
            case it is resolved using online services or as the appropriate astropy SkyCoord object.
            ICRS coordinates may also be entered as a string.
        radius : `number`
            The radius of the circular region to query.

        Returns
        -------
        catalogue : `astropy.table.Table`
            An astropy Table with all the information about the star

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
        """

        Parameters
        ----------
        coord : `str`, `astropy.coordinates.SkyCoord`
            The target around which to search. It may be specified as a string in which
            case it is resolved using online services or as the appropriate astropy SkyCoord object.
            ICRS coordinates may also be entered as a string.
        radius : `number`
            The radius of the circular region to query.
        width : `number`
            The width in  of the square region to query.
        height : `number`
            When set in addition to ``width``, the queried region becomes
            rectangular, with the specified ``width`` and ``height``.
        columns : `list`
            List of strings with the keyword to fetch the catalog.
            If ``columns=None`` it will download all the columns.
            If ``columns="simple"`` it will download only the columns for
            the code, epoch and astrometric parameters.
        row_limit : int
            Maximum number of rows that will be fetched from the result
            (set to -1 for unlimited). Default: ``row_limit=10_000_000``
        timeout : `number`
            timeout for connecting to server in seconds. Default: ``timeout=600``
        **kwargs
            Any other keyword argument will be passed to astroquery.vizier.Vizier

        Returns
        -------
        catalogue : `astropy.table.Table`
            An astropy Table with all the information about the star
        """
        if columns == 'simple':
            columns = self.get_simple_columns()
        elif columns is None:
            columns = ['**']
        vquery = Vizier(columns=columns, row_limit=row_limit, timeout=timeout, **kwargs)
        catalogue = vquery.query_region(coord, radius=radius, width=width, height=height, catalog=self.cat_path, cache=False)
        return catalogue

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        return f'<VizierCatalogue: {self.name} defined in https://vizier.cds.unistra.fr/viz-bin/VizieR-3?-source={self.cat_path}>'


class LineaGaiaCatalogue(Catalogue):
    """Gaia catalogue served by LIneA TAP."""

    def __init__(self, tap_url='https://userquery.linea.org.br/tap', language='PostgreSQL', **kwargs):
        self.tap_url = tap_url
        self.language = language
        super(LineaGaiaCatalogue, self).__init__(**kwargs)

    @staticmethod
    def _format_value(value):
        text = str(value).strip()
        if text.isdigit():
            return text
        return "'{}'".format(text.replace("'", "''"))

    def _run_query(self, query):
        session = requests.Session()
        tap = pyvo.dal.TAPService(self.tap_url, session=session)
        return tap.run_sync(query, language=self.language).to_table()

    def _format_columns(self, columns):
        if columns == 'simple':
            columns = self.get_simple_columns()
        if columns is None:
            return '*'
        return ', '.join(dict.fromkeys(columns))

    def _format_column_filters(self, column_filters):
        if not column_filters:
            return []
        filters = []
        for column, expression in column_filters.items():
            filters.append(f'{column} {str(expression).strip()}')
        return filters

    def search_star(self, code=None, coord=None, radius=None):
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
        del timeout, kwargs

        if not isinstance(coord, SkyCoord):
            coord = SkyCoord(coord, unit=('hourangle', 'deg'))
        coord = coord.icrs

        if radius is not None:
            radius = u.Quantity(radius).to(u.deg)
            region_filter = f'q3c_radial_query({self.ra}, {self.dec}, {coord.ra.deg}, {coord.dec.deg}, {radius.value})'
        elif width is not None:
            width = u.Quantity(width).to(u.deg)
            height = u.Quantity(height if height is not None else width).to(u.deg)
            half_width = width.value / 2.0
            half_height = height.value / 2.0
            ra0 = coord.ra.wrap_at(360 * u.deg).deg
            dec0 = coord.dec.deg
            ra_min = (ra0 - half_width) % 360.0
            ra_max = (ra0 + half_width) % 360.0
            vertices = [
                ra_min, dec0 - half_height,
                ra_max, dec0 - half_height,
                ra_max, dec0 + half_height,
                ra_min, dec0 + half_height,
            ]
            vertex_string = ', '.join(f'{value:.16f}' for value in vertices)
            region_filter = (
                f'q3c_poly_query({self.ra}, {self.dec}, ARRAY[{vertex_string}]::double precision[])'
            )
        else:
            raise ValueError('At least a radius or width should be given as input')

        query_filters = [region_filter]
        query_filters.extend(self._format_column_filters(column_filters))
        limit = f' LIMIT {int(row_limit)}' if row_limit and row_limit > 0 else ''
        query = f"SELECT {self._format_columns(columns)} FROM {self.cat_path} WHERE {' AND '.join(query_filters)}{limit}"
        catalogue = self._run_query(query)
        return [] if len(catalogue) == 0 else [catalogue]

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        return f'<LineaGaiaCatalogue: {self.name} defined in {self.tap_url} ({self.cat_path})>'


def is_timeout_error(exc):
    """Returns True when an exception indicates a timeout condition."""
    if isinstance(exc, (requests.exceptions.Timeout, TimeoutError)):
        return True
    text = str(exc).lower()
    return 'timed out' in text or 'timeout' in text


def should_fallback_to_gaiadr3(catalogue, exc):
    """Returns True when a TapLinea query should fallback to VizieR Gaia DR3."""
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

taplinea = LineaGaiaCatalogue(name='GaiaDR3', cat_path='gaia_dr3.source', code='source_id', ra='ra', dec='dec',
                              pmra='pmra', pmdec='pmdec', epoch='ref_epoch', parallax='parallax',
                              rad_vel='radial_velocity',
                              band={'G': 'phot_g_mean_mag', 'BP': 'phot_bp_mean_mag', 'RP': 'phot_rp_mean_mag'},
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
lineagaia = taplinea

epoch_ubsc = Time('J1991.25', scale='tdb')
ubsc = VizierCatalogue(name='UBSC', cat_path='J/AJ/164/36/table2', code='HIP', ra='RA_ICRS', dec='DE_ICRS',
                       pmra='pmRA', pmdec='pmDE', epoch=epoch_ubsc, parallax='plx', band={'Hp': 'Hpmag'},
                       errors=['e_RA_ICRS', 'e_DE_ICRS', 'e_pmRA', 'e_pmDE', 'e_plx', None])

allowed_catalogues = SelectDefault(instance=Catalogue,
                                   defaults={'TapLinea': taplinea, 'taplinea': taplinea, 'lineagaia': taplinea,
                                             'lineagaiadr3': taplinea, 'gaiadr2': gaiadr2, 'gaiaedr3': gaiaedr3,
                                             'gaiadr3': gaiadr3, 'ubsc': ubsc})
