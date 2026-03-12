import astropy.units as u
import numpy as np

__all__ = ['search_code_mpc']

class MPCServiceError(ConnectionError):
    """Raised when the MPC observatory service cannot be reached."""


def search_code_mpc(code, timeout=30):
    """Reads one observatory from the MPC database using the MPC API.

    Parameters
    ----------
    code : `str`
        MPC observatory code.
    timeout : `int`, optional
        Timeout for the HTTP request, in seconds.

    Returns
    -------
    name, site : `tuple`
        Observatory name and site as an Astropy EarthLocation object.

    Raises
    ------
    ValueError
        If the observatory code is not found in the MPC database.
    MPCServiceError
        If the MPC service cannot be reached.
    """
    import requests
    from astropy.coordinates import EarthLocation

    code = str(code).strip()

    url = "https://data.minorplanetcenter.net/api/obscodes"
    payload = {"obscode": code, "format": "JSON"}
    headers = {"User-Agent": "SORA observatory lookup"}

    try:
        response = requests.get(url, json=payload, headers=headers, timeout=timeout)
        response.raise_for_status()
        data = response.json()
    except requests.RequestException as e:
        raise MPCServiceError(
            f"Could not access MPC Observatory Codes API: {e}"
        ) from e

    if not data:
        raise ValueError(f"code {code} could not be located in MPC database")

    name = data.get("name_utf8") or data.get("name") or code

    lon = float(data["longitude"]) * u.deg
    rcphi = float(data["rhocosphi"]) * 6378.137 * u.km
    rsphi = float(data["rhosinphi"]) * 6378.137 * u.km

    site = EarthLocation.from_geocentric(
        rcphi * np.cos(lon),
        rcphi * np.sin(lon),
        rsphi,
    )

    return name, site