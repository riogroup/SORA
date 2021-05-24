import astropy.units as u
import numpy as np

__all__ = ['search_code_mpc']


def search_code_mpc():
    """Reads the MPC Observer Database.

    Returns
    -------
    observatories : `dict`
        A python dictionary with all the sites as Astropy EarthLocation objects.
    """
    from astropy.coordinates import EarthLocation
    from astroquery.mpc import MPC

    obs = MPC.get_observatory_codes()
    observatories = {}
    for line in obs:
        code = line['Code']
        lon = line['Longitude'] * u.deg
        rcphi = line['cos'] * 6378.137 * u.km
        rsphi = line['sin'] * 6378.137 * u.km
        name = line['Name']
        site = EarthLocation.from_geocentric(rcphi * np.cos(lon), rcphi * np.sin(lon), rsphi)
        observatories[code] = (name, site)
    return observatories
