"""Utility functions for file, Horizons, and SPICE-kernel ephemerides."""

import warnings

import requests
import pathlib
import base64
import http
import re
from datetime import datetime as dt

satellites_bsp = {
    "PLUTO": "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/satellites/plu060.bsp",
    "1930 BM": "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/satellites/plu060.bsp",
}


def _resolve_horizons_options(timeout=None, cache=None):
    """Resolve explicit Horizons options against configured defaults."""
    if timeout is None or cache is None:
        from sora.config import get_config

        services = get_config().services
        if timeout is None:
            timeout = services.horizons_timeout
        if cache is None:
            cache = services.horizons_cache
    return timeout, cache


def getBSPfromJPL(identifier, initial_date, final_date, email=None,
                  directory='./', timeout=None):
    """
    Retrieve SPK (Spacecraft and Planet Kernel) files from JPL Horizons API for given object identifiers.

    Parameters
    ----------
    identifier : str or list of str
        A single object identifier or a list of identifiers: names, numbers, provisional designations, and SPK ID.
        If an identifier is an SPKID it MUST start with the prefix 'SPK' (case-insensitive), and followed by numeric digits;
        the function strips the 'SPK' prefix and any non-digit characters, using only the digits for the lookup.
        Other identifiers are used as provided.
    initial_date : str
        Start date in 'YYYY-MM-DD' format. Must be at least 32 days before `final_date`.
    final_date : str
        End date in 'YYYY-MM-DD' format. Must be at least 32 days after `initial_date`.
    email : str, optional
        User email for JPL API (reserved for backward compatibility; not used).
    directory : str, optional
        Path to the directory where BSP files will be saved. Default is current directory.
    timeout : int, float, optional
        HTTP timeout in seconds. When omitted, uses the configured Horizons
        timeout.

    Returns
    -------
    None
        Writes BSP files to the specified directory.

    Raises
    ------
    ValueError
        If date range is invalid, directory does not exist, or an SPK identifier is malformed (missing digits after 'SPK').
    Exception
        For HTTP errors, missing SPK generation, or file writing issues.
    """
    timeout, _ = _resolve_horizons_options(timeout=timeout, cache=False)
    ids = [identifier] if isinstance(identifier, str) else list(identifier)
    ids = [str(obj) for obj in ids]

    date1 = dt.strptime(initial_date, "%Y-%m-%d")
    date2 = dt.strptime(final_date,  "%Y-%m-%d")
    if (date2 - date1).days <= 32:
        raise ValueError("final_date must be more than 32 days after initial_date")

    path = pathlib.Path(directory)
    if not path.exists():
        raise ValueError(f"Directory {path} does not exist")

    base_url = "https://ssd.jpl.nasa.gov/api/horizons.api"
    success, failures = 0, []
    max_len = max(len(obj) for obj in ids)

    def _fetch_data(params):
        """
        Fetch data from the JPL Horizons API with given parameters.

        Parameters
        ----------
        params : dict
            Dictionary of query parameters for the API request.

        Returns
        -------
        requests.Response
            The HTTP response object from the GET request.

        Raises
        ------
        Exception
            If the HTTP status code is not 200 (OK).
        """
        r = requests.get(
            base_url,
            params=params,
            stream=True,
            timeout=timeout,
        )
        if r.status_code != 200:
            raise Exception(f"HTTP {r.status_code} - {http.HTTPStatus(r.status_code).phrase}")
        return r

    for obj in ids:
        # If identifier starts with 'SPK', extract only digits after it
        spk = False
        if obj.upper().startswith('SPK'):
            digits = re.findall(r"\d+", obj)
            if not digits:
                raise ValueError(f"No digits found in SPK identifier '{obj}'")
            obj = ''.join(digits)
            fname   = f"SPKID{obj.replace(' ', '')}.bsp"
            spk = True
        else:
            fname   = f"{obj.replace(' ', '')}.bsp"

        status = f"Retrieving {obj:{max_len}} …"
        print(status, end="\r")

        params = {
            "format":     "json",
            "EPHEM_TYPE": "SPK",
            "OBJ_DATA":   "YES",
            "COMMAND":    f"'DES={obj}'",
            "START_TIME": date1.strftime("%Y-%b-%d"),
            "STOP_TIME":  date2.strftime("%Y-%b-%d"),
        }

        try:
            r = _fetch_data(params)
            result = r.json().get("result", "")

            if ("No matches found." in result or "Comet AND asteroid index search" in result) and not spk:
                params["COMMAND"] = f"'{obj};'"
                r = _fetch_data(params)

            data = r.json()
            if "error" in data and "SPK creation is not available" in data["error"]:
                url = satellites_bsp.get(obj.upper())
                if not url:
                    raise Exception("No fallback URL for satellite")
                r = requests.get(url, stream=True, timeout=timeout)
                r.raise_for_status()
                content = r.content
            elif "spk" in data:
                content = base64.b64decode(data["spk"])
            else:
                raise Exception("SPK file not generated")


            outpath = path / fname
            with open(outpath, "wb") as f:
                f.write(content)
            success += 1

        except Exception as e:
            failures.append(f"{obj}: {e}")

        finally:
            print(" " * (len(status) + 5), end="\r")

    summary = f"Done — successes: {success}, failures: {len(failures)}"
    print(summary.ljust(len(summary) + max_len))
    if failures:
        print("Failures detail:")
        for line in failures:
            print("  — ", line)
        print("     Check if the identifier(s) is(are) correct.")



def ephem_kernel(time, target, observer, kernels, output='ephemeris'):
    """Calculates the ephemeris from kernel files.

    Parameters
    ----------
    time : `str`, `astropy.time.Time`
        Reference instant to calculate ephemeris. It can be a string
        in the ISO format (yyyy-mm-dd hh:mm:ss.s) or an astropy Time object.

    target : `str`
        IAU (kernel) code of the target.

    observer : `str`, `sora.Observer`, `sora.Spacecraft`
        IAU (kernel) code of the observer, a SORA observer object, a SORA
        spacecraft object, or one of ``'geocenter'`` and ``'barycenter'``.
        String IAU codes must be present in the loaded kernels.

    kernels : `list`, `str`
        List of paths for all the kernels.

    output : `str`, optional, default='ephemeris'
        The output of data. ``ephemeris`` will output the observed position,
        while ``vector`` will output the Cartesian state vector, without
        light time correction.

    Returns
    -------
    coord : `astropy.coordinates.SkyCoord`
        ICRS coordinate of the target when ``output='ephemeris'``. Cartesian
        state vector of the target relative to the observer when
        ``output='vector'``.
    """
    import numpy as np
    import astropy.units as u
    import astropy.constants as const
    import spiceypy as spice

    from astropy.coordinates import SkyCoord, GCRS
    from astropy.time import Time
    from astropy.utils.exceptions import AstropyWarning
    from sora.observer import Observer, Spacecraft

    origins = {'geocenter': '399', 'barycenter': '0'}
    location = origins.get(observer)
    if not location and isinstance(observer, str):
        location = observer
    if isinstance(observer, (Observer, Spacecraft)):
        location = str(getattr(observer, "spkid", None))
    if not location:
        raise ValueError("observer must be 'geocenter', 'barycenter' or an observer object.")
    if output not in ['ephemeris', 'vector']:
        raise ValueError("output must be 'ephemeris' or 'vector'")

    if type(kernels) == str:
        kernels = [kernels]
    for kern in kernels:
        spice.furnsh(kern)
    time = Time(time)
    t0 = Time('J2000', scale='tdb')
    dt = (time - t0)
    delt = 0 * u.s

    # calculates vector Solar System Barycenter -> Observer
    if isinstance(observer, Observer):
        position1 = spice.spkpos(location, dt.sec, 'J2000', 'NONE', '0')[0]
        position1 = SkyCoord(*position1.T * u.km, representation_type='cartesian')
        itrs = observer.site.get_itrs(obstime=time)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', AstropyWarning)
            gcrs = itrs.transform_to(GCRS(obstime=time))
        position1 = SkyCoord(position1.cartesian + gcrs.cartesian, representation_type='cartesian')
    elif isinstance(observer, Spacecraft):
        spice.kclear()  # necessary because observer.get_vector() may load different kernels
        position1 = observer.get_vector(time=time, origin='barycenter')
        for kern in kernels:
            spice.furnsh(kern)
    else:
        position1 = spice.spkpos(location, dt.sec, 'J2000', 'NONE', '0')[0]
        position1 = SkyCoord(*position1.T * u.km, representation_type='cartesian')

    while True:
        # calculates new time
        tempo = dt - delt
        # calculates vector Solar System Barycenter -> Object
        position2 = spice.spkpos(target, tempo.sec, 'J2000', 'NONE', '0')[0]
        position2 = SkyCoord(*position2.T * u.km, representation_type='cartesian')
        position = position2.cartesian - position1.cartesian
        # calculates new light time
        delt = (position.norm() / const.c).decompose()
        # if difference between new and previous light time is smaller than 0.001 sec, then continue.
        if output == 'vector' or np.all(np.absolute(((dt - tempo) - delt).sec) < 0.001):
            break
    coord = SkyCoord(position, representation_type='cartesian')
    spice.kclear()
    if output == 'ephemeris':
        coord = SkyCoord(ra=coord.spherical.lon, dec=coord.spherical.lat,
                         distance=coord.spherical.distance, obstime=time)
    if not coord.isscalar and len(coord) == 1:
        coord = coord[0]
    return coord


def ephem_horizons(time, target, observer, id_type='smallbody',
                   output='ephemeris', timeout=None, cache=None):
    """Calculates the ephemeris from Horizons.

    Parameters
    ----------
    time : `str`, `astropy.time.Time`
        Reference instant to calculate ephemeris. It can be a string
        in the ISO format (yyyy-mm-dd hh:mm:ss.s) or an astropy Time object.

    target : `str`
        Target name or identifier accepted by Horizons.

    observer : `str`, `sora.Observer`, `sora.Spacecraft`
        Horizons observer code, a SORA observer object, a SORA spacecraft
        object, or one of ``'geocenter'`` and ``'barycenter'``.

    id_type : `str`, optional, default='smallbody'
        Type of target object options: ``smallbody``, ``majorbody`` (planets but
        also anything that is not a small body), ``designation``, ``name``,
        ``asteroid_name``, ``comet_name``, or ``id`` (Horizons id number).

    output : `str`, optional, default='ephemeris'
        The output of data. ``ephemeris`` will output the observed position,
        while ``vector`` will output the Cartesian state vector, without
        light time correction.

    timeout : `int`, `float`, optional
        Horizons connection timeout in seconds. When omitted, uses the
        configured value.

    cache : `bool`, optional
        Whether Astroquery may reuse cached Horizons responses. When omitted,
        uses the configured value.

    Returns
    -------
    coord : `astropy.coordinates.SkyCoord`
        ICRS coordinate of the target when ``output='ephemeris'``. Cartesian
        state vector from Horizons when ``output='vector'``.

    Notes
    -----
    If the interval of time is larger than 30 days or so, a timeout error may be raised.
    The maximum interval will depend on the user connection.
    """
    import astropy.units as u

    from astroquery.jplhorizons import Horizons
    from astropy.time import Time
    from astropy.coordinates import SkyCoord
    from sora.observer import Observer, Spacecraft
    from scipy import interpolate

    origins = {'geocenter': '@399', 'barycenter': '@0'}
    location = origins.get(observer)
    if not location and isinstance(observer, str):
        location = observer
    if isinstance(observer, Observer):
        if getattr(observer, "code", None) is None:
            location = {'lon': observer.lon.deg, 'lat': observer.lat.deg, 'elevation': observer.height.to(u.km).value}
        else:
            location = f'{getattr(observer, "code", "")}@{getattr(observer, "spkid", "")}'
    elif isinstance(observer, Spacecraft):
        location = f"500@{getattr(observer, 'spkid', '')}"
    if not location:
        raise ValueError("observer must be 'geocenter', 'barycenter' or an observer object.")
    if output not in ['ephemeris', 'vector']:
        raise ValueError("output must be 'ephemeris' or 'vector'")

    timeout, cache = _resolve_horizons_options(
        timeout=timeout,
        cache=cache,
    )
    time = Time(time)
    time1 = getattr(time, {'ephemeris': 'utc', 'vector': 'tdb'}[output]).jd
    if not time.isscalar and len(time) > 50:
        step = '10m'
        if time.max() - time.min() > 30 * u.day:
            warnings.warn('Time interval may be too long. A timeout error may be raised.')
        if time.max() - time.min() <= 1 * u.day:
            step = '1m'
        time2 = {'start': (time.min() - 10 * u.min).iso.split('.')[0],
                 'stop': (time.max() + 10 * u.min).iso.split('.')[0],
                 'step': step,
                 }
    else:
        time2 = time1

    if getattr(observer, 'ephem', None) not in ['horizons', None]:
        warnings.warn('Ephemeris using kernel for the observer and Horizons for the target is under construction. '
                      'We will use only Horizons.')
    id_type = None if id_type == 'majorbody' else id_type
    ob = Horizons(id=target, id_type=id_type, location=location, epochs=time2)
    ob.TIMEOUT = timeout

    if output == 'ephemeris':
        eph = ob.ephemerides(extra_precision=True, cache=cache)
        obstime = Time(eph['datetime_jd'], format='jd', scale='utc')
        pos = SkyCoord(eph['RA'], eph['DEC'], eph['delta'], frame='icrs', obstime=obstime)
    else:
        vec = ob.vectors(refplane='earth', cache=cache)
        obstime = Time(vec['datetime_jd'], format='jd', scale='tdb')
        pos = SkyCoord(*[vec[i] for i in ['x', 'y', 'z']] * u.AU, representation_type='cartesian', obstime=obstime)

    if isinstance(time2, dict):
        spl_x = interpolate.CubicSpline(obstime.jd, pos.cartesian.x.to(u.km))
        spl_y = interpolate.CubicSpline(obstime.jd, pos.cartesian.y.to(u.km))
        spl_z = interpolate.CubicSpline(obstime.jd, pos.cartesian.z.to(u.km))
        pos = SkyCoord(x=spl_x(time1), y=spl_y(time1), z=spl_z(time1), unit=u.km, representation_type='cartesian')

    if output == 'ephemeris':
        pos = SkyCoord(ra=pos.spherical.lon, dec=pos.spherical.lat, distance=pos.spherical.distance)

    if not pos.isscalar and len(pos) == 1:
        pos = pos[0]

    return pos
