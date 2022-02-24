from abc import ABCMeta, abstractmethod
from dataclasses import dataclass

import astropy.units as u
import numpy as np
from astropy.time import Time
from astroquery.vizier import Vizier


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

    @abstractmethod
    def search_star(self):
        pass

    @abstractmethod
    def search_region(self, coord, radius=None, width=None, height=None, columns=None, verbose=False,
                      row_limit=10000000, timeout=600, **kwargs):
        pass

    def parse_catalogue(self, table):
        """

        Parameters
        ----------
        table : `astropy.table.Table`

        Returns
        -------

        """
        output: dict = {'code': table[self.code].tolist()}
        params = ['ra', 'dec', 'pmra', 'pmdec', 'parallax', 'rad_vel']
        units = [u.deg, u.deg, u.mas / u.year, u.mas / u.year, u.mas, u.km / u.s]
        for p, unit in zip(params, units):
            if not getattr(self, p, None):
                vals = np.zeros(len(table)) * unit
            else:
                vals = table[getattr(self, p).replace('(', '_').replace(')', '_')].quantity
            vals[np.where(np.isnan(vals))] = 0 * unit
            output[p] = vals
        if isinstance(self.epoch, Time):
            output['epoch'] = Time(np.repeat(self.epoch, len(table)))
        else:
            output['epoch'] = Time(table[self.epoch].quantity, format='jyear')
        if self.band is not None:
            bands = {keys: table[item].quantity for keys, item in self.band.items() if item in table.columns}
            output['band'] = bands
        return output


class VizierCatalogue(Catalogue):

    def search_star(self, code=None, coord=None, radius=None, verbose=False):
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

    def search_region(self, coord, radius=None, width=None, height=None, columns=None, verbose=False,
                      row_limit=10000000, timeout=600, **kwargs):
        col_keys = ['code', 'ra', 'dec', 'pmra', 'pmdec', 'parallax', 'rad_vel', 'epoch']
        if columns == 'simple':
            columns = [getattr(self, p) for p in col_keys if hasattr(self, p) and isinstance(getattr(self, p), str)] + \
                      list(self.band.values() if hasattr(self, 'band') else {})
        elif columns is None:
            columns = ['**']
        vquery = Vizier(columns=columns, row_limit=row_limit, timeout=timeout, **kwargs)
        catalogue = vquery.query_region(coord, radius=radius, width=width, height=height, catalog=self.cat_path, cache=False)
        return catalogue


gaiadr2 = VizierCatalogue(name='GaiaDR2', cat_path='I/345/gaia2', code='Source', ra='RA_ICRS', dec='DE_ICRS',
                          pmra='pmRA', pmdec='pmDE', epoch='Epoch', parallax='Plx', rad_vel='RV', band={'G': 'Gmag'})

gaiaedr3 = VizierCatalogue(name='GaiaEDR3', cat_path='I/350/gaiaedr3', code='Source', ra='RA_ICRS', dec='DE_ICRS',
                           pmra='pmRA', pmdec='pmDE', epoch='Epoch', parallax='Plx', rad_vel='RVDR2', band={'G': 'Gmag'})


catalogs = {'gaiadr2': gaiadr2, 'gaiaedr3': gaiaedr3}
