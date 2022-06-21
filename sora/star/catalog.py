from abc import ABCMeta, abstractmethod
from dataclasses import dataclass

import astropy.units as u
import numpy as np
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
                vals = table[getattr(self, p).replace('(', '_').replace(')', '_')].quantity
            vals[np.where(np.isnan(vals))] = 0 * unit
            output[p] = vals
            if err is not None and err in table.colnames:
                vals = table[err].quantity
                vals[np.where(np.isnan(vals))] = 0 * unit
                errors[p] = vals
        output['errors'] = errors
        if isinstance(self.epoch, Time):
            output['epoch'] = Time(np.repeat(self.epoch, len(table)))
        else:
            output['epoch'] = Time(table[self.epoch].quantity, format='jyear')
        if self.band is not None:
            bands = {keys: table[item].quantity for keys, item in self.band.items() if item in table.columns}
            output['band'] = bands
        return output


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
        col_keys = ['code', 'ra', 'dec', 'pmra', 'pmdec', 'parallax', 'rad_vel', 'epoch']
        if columns == 'simple':
            columns = [getattr(self, p) for p in col_keys if hasattr(self, p) and isinstance(getattr(self, p), str)] + \
                      list(self.band.values() if hasattr(self, 'band') else {})
        elif columns is None:
            columns = ['**']
        vquery = Vizier(columns=columns, row_limit=row_limit, timeout=timeout, **kwargs)
        catalogue = vquery.query_region(coord, radius=radius, width=width, height=height, catalog=self.cat_path, cache=False)
        return catalogue

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        return f'<VizierCatalogue: {self.name} defined in https://vizier.cds.unistra.fr/viz-bin/VizieR-3?-source={self.cat_path}>'


gaiadr2 = VizierCatalogue(name='GaiaDR2', cat_path='I/345/gaia2', code='Source', ra='RA_ICRS', dec='DE_ICRS',
                          pmra='pmRA', pmdec='pmDE', epoch='Epoch', parallax='Plx', rad_vel='RV', band={'G': 'Gmag'},
                          errors=['e_RA_ICRS', 'e_DE_ICRS', 'e_pmRA', 'e_pmDE', 'e_Plx', 'e_RV'])

gaiaedr3 = VizierCatalogue(name='GaiaEDR3', cat_path='I/350/gaiaedr3', code='Source', ra='RA_ICRS', dec='DE_ICRS',
                           pmra='pmRA', pmdec='pmDE', epoch='Epoch', parallax='Plx', rad_vel='RVDR2', band={'G': 'Gmag'},
                           errors=['e_RA_ICRS', 'e_DE_ICRS', 'e_pmRA', 'e_pmDE', 'e_Plx', 'e_RVDR2'])

epoch = Time('J2016.0', scale='tdb')
gaiadr3 = VizierCatalogue(name='GaiaDR3', cat_path='I/355/gaiadr3', code='Source', ra='RA_ICRS', dec='DE_ICRS',
                          pmra='pmRA', pmdec='pmDE', epoch=epoch, parallax='Plx', rad_vel='RV', band={'G': 'Gmag'},
                          errors=['e_RA_ICRS', 'e_DE_ICRS', 'e_pmRA', 'e_pmDE', 'e_Plx', 'e_RV'])


allowed_catalogues = SelectDefault(instance=VizierCatalogue,
                                   defaults={'gaiadr2': gaiadr2, 'gaiaedr3': gaiaedr3, 'gaiadr3': gaiadr3})
