from os import PathLike, fspath
from pathlib import Path

import astropy.units as u
import numpy as np
import pyvo
import requests

from sora.config import get_config

__all__ = ['search_code_mpc']

_MPC_CODE_CACHE = {}


def _resolve_observer_ephem(ephem):
    """Resolve observer ephemerides to Horizons or an ordered kernel list.

    Bare names are looked up in the planetary catalogue, including locally
    registered kernels. Paths are passed through to SPICE. Multipart entries
    expand in place, preserving their order relative to other kernels.
    """
    if isinstance(ephem, str) and ephem.casefold() == 'horizons':
        return 'horizons'
    if isinstance(ephem, (str, PathLike)):
        ephem = [ephem]
    if not isinstance(ephem, list):
        raise TypeError('ephem must be "horizons", a kernel name or path, or a list of kernels')
    if not ephem:
        raise ValueError('ephem must contain at least one kernel')
    if any(
        not isinstance(kernel, (str, PathLike))
        or not isinstance(fspath(kernel), str)
        or not fspath(kernel).strip()
        for kernel in ephem
    ):
        raise TypeError('Kernel names and paths must be non-empty strings or paths')

    database = None
    kernels = []
    for kernel in ephem:
        path = fspath(kernel)
        # Explicit paths do not need a catalogue lookup or network access.
        if (
            isinstance(kernel, PathLike)
            or '/' in path
            or '\\' in path
            or Path(path).suffix
            or Path(path).is_file()
        ):
            kernels.append(path)
            continue
        if database is None:
            from sora.ephem.planetary_kernels import PlanetaryKernelDB

            database = PlanetaryKernelDB(config=get_config())
        kernels.extend(database.get_planetary_kernels(path))
    return kernels


def search_code_mpc(code):
    """Queries the Minor Planet Center (MPC) Observer codes SBN mirror.

    Parameters
    ----------
    code : `str`, `int`
        MPC observatory code.

    Returns
    -------
    name : `str`
        Observatory name from the MPC database.

    site : `astropy.coordinates.EarthLocation`
        Observatory site as an Astropy EarthLocation object.

    Raises
    ------
    ValueError
        Raised when the MPC code is not found in the database.

    Notes
    -----
    Query results are cached in memory by MPC code.
    """
    from astropy.coordinates import EarthLocation
    import warnings

    code = str(code).strip()
    if code in _MPC_CODE_CACHE:
        return _MPC_CODE_CACHE[code]
    
    url = get_config().services.linea_tap_url
    session = requests.Session()
    tap = pyvo.dal.TAPService(url, session=session)
    
    query = f"SELECT * FROM mpc_sbn.obscodes WHERE obscode = '{code}'"
    warnings.warn(f'Querying code {code} in the Linea MPC Observer Database...')
    result = tap.run_sync(query)
    table = result.to_table()

    if len(table) == 0:
        raise ValueError(f'code {code} could not be located in MPC database')

    line = table[0]
    lon = line['longitude'] * u.deg
    rcphi = line['rhocosphi'] * 6378.137 * u.km
    rsphi = line['rhosinphi'] * 6378.137 * u.km
    name = line['name']
    site = EarthLocation.from_geocentric(rcphi * np.cos(lon),
                                         rcphi * np.sin(lon),
                                         rsphi)
    _MPC_CODE_CACHE[code] = (name, site)
    return _MPC_CODE_CACHE[code]
