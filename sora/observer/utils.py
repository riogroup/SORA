import astropy.units as u
import numpy as np
import pyvo
import requests

__all__ = ['search_code_mpc']

_MPC_CODE_CACHE = {}


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
    
    url = "https://userquery.linea.org.br/tap"
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
