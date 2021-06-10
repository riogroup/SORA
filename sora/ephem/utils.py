import warnings

__all__ = ['getBSPfromJPL', 'ephem_kernel', 'ephem_horizons']


def getBSPfromJPL(identifier, initial_date, final_date, email, directory='./'):
    """Downloads BSP files from JPL database.

    BSP files, which have information to generate the ephemeris of the objects,
    will be downloaded and named as (without spaces): '[identifier].bsp'.

    Important
    ---------
    It is not possible to download BSP files for Planets or Satellites.

    Parameters
    ----------
    identifier : `str`, `list`
        Identifier of the object(s).
        It can be the `name`, `number` or `SPK ID`.
        It can also be a list of objects.

        Examples: ``'2137295'``, ``'1999 RB216'``, ``'137295'``,
        ``['2137295', '136199', '1999 RC216', 'Chariklo']``.

    initial_date : `str`
        Date the bsp file is to begin, within span `1900-2100`.

        Examples: ``'2003-02-01'``, ``'2003-3-5'``.

    final_date : `str`
        Date the bsp file is to end, within span [1900-2100].
        Must be more than 32 days later than `initial_date`.

        Examples: ``'2006-01-12'``, ``'2006-1-12'``.

    email : `str`
        User's e-mail contact address.
        Required by JPL web service.

        Example: ``username@user.domain.name``.

    directory : `str`
        Directory path to save the bsp files.
    """
    import pathlib
    import shutil
    import requests
    from datetime import datetime

    date1 = datetime.strptime(initial_date, '%Y-%m-%d')
    date2 = datetime.strptime(final_date, '%Y-%m-%d')
    diff = date2 - date1

    if diff.days <= 32:
        raise ValueError('The [final_date] must be more than 32 days later than [initial_date]')
    else:
        if isinstance(identifier, str):
            identifier = [identifier]

        url_jpl = 'https://ssd.jpl.nasa.gov/x/smb_spk.cgi?OPTION=Make+SPK'
        lim, opt, n = 10, 'y', len(identifier)

        if n > lim:
            parameters = {'OBJECT': identifier[0], 'START': date1.strftime('%Y-%b-%d'),
                          'STOP': date2.strftime('%Y-%b-%d'), 'EMAIL': email, 'TYPE': '-B'}

            t0 = datetime.now()
            r = requests.get(url_jpl, params=parameters, stream=True)
            tf = datetime.now()

            bsp_format = r.headers['Content-Type']
            if r.status_code == requests.codes.ok and bsp_format == 'application/download':
                size0 = int(r.headers["Content-Length"]) / 1024 / 1024  # MB

                print('Estimated time to download {} ({:.3f} MB) files: {}'.
                      format(n, n * size0, n * (tf - t0)))

                opt = input('\nAre you sure? (y/n):')
            else:
                raise ValueError('It was not able to download the bsp file for object {}\n'.
                                 format(identifier[0]))

        if opt in ['y', 'Y', 'YES', 'yes']:
            if directory != './':
                path = pathlib.Path(directory)
                if not path.exists():
                    raise ValueError('The directory {} does not exist!'.format(path))
                directory += '/'

            print("\nDownloading bsp file(s) ...\n")

            m, size = 0, 0.0
            failed = []

            t0 = datetime.now()
            for obj in identifier:
                filename = obj.replace(' ', '') + '.bsp'

                parameters = {'OBJECT': obj, 'START': date1.strftime('%Y-%b-%d'),
                              'STOP': date2.strftime('%Y-%b-%d'), 'EMAIL': email, 'TYPE': '-B'}

                r = requests.get(url_jpl, params=parameters, stream=True)
                bsp_format = r.headers['Content-Type']
                if r.status_code == requests.codes.ok and bsp_format == 'application/download':
                    size += int(r.headers["Content-Length"]) / 1024 / 1024
                    m += 1
                    with open(directory + filename, 'wb') as f:
                        r.raw.decode_content = True
                        shutil.copyfileobj(r.raw, f)
                else:
                    failed.append(obj)
            tf = datetime.now()

            print("{} ({:.3f} MB) file(s) was/were downloaded".format(m, size))
            print("Download time: {}".format(tf - t0))

            if len(failed) > 0:
                raise ValueError(
                    'It was not able to download the bsp files for next objects: {}'.format(sorted(failed)))


def ephem_kernel(time, target, observer, kernels, output='ephemeris'):
    """Calculates the ephemeris from kernel files.

    Parameters
    ----------
    time : `str`, `astropy.time.Time`
        Reference instant to calculate ephemeris. It can be a string
        in the ISO format (yyyy-mm-dd hh:mm:ss.s) or an astropy Time object.

    target : `str`
        IAU (kernel) code of the target.

    observer : `str`
        IAU (kernel) code of the observer.

    kernels : `list`, `str`
        List of paths for all the kernels.

    output : `str`
        The output of data. ``ephemeris`` will output the observed position,
        while ``vector`` will output the Cartesian state vector, without
        light time correction.

    Returns
    -------
    coord : `astropy.coordinates.SkyCoord`
        ICRS coordinate of the target.
    """
    import numpy as np
    import astropy.units as u
    import astropy.constants as const
    import spiceypy as spice

    from astropy.coordinates import SkyCoord
    from astropy.time import Time
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
    if isinstance(observer, (Observer, Spacecraft)):
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


def ephem_horizons(time, target, observer, id_type='smallbody', output='ephemeris'):
    """Calculates the ephemeris from Horizons.

    Parameters
    ----------
    time : `str`, `astropy.time.Time`
        Reference instant to calculate ephemeris. It can be a string
        in the ISO format (yyyy-mm-dd hh:mm:ss.s) or an astropy Time object.

    target : `str`
        IAU (kernel) code of the target.

    observer : `str`
        IAU (kernel) code of the observer.

    id_type : `str`
        Type of target object options: ``smallbody``, ``majorbody`` (planets but
        also anything that is not a small body), ``designation``, ``name``,
        ``asteroid_name``, ``comet_name``, ``id`` (Horizons id number), or
        ``smallbody`` (find the closest match under any id_type).

    output : `str`
        The output of data. ``ephemeris`` will output the observed position,
        while ``vector`` will output the Cartesian state vector, without
        light time correction.

    Returns
    -------
    coord : `astropy.coordinates.SkyCoord`
        ICRS coordinate of the target.

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
    if isinstance(observer, (Observer, Spacecraft)):
        location = f'{getattr(observer, "code", "")}@{getattr(observer, "spkid", "")}'
    if not location:
        raise ValueError("observer must be 'geocenter', 'barycenter' or an observer object.")
    if output not in ['ephemeris', 'vector']:
        raise ValueError("output must be 'ephemeris' or 'vector'")

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
    ob = Horizons(id=target, id_type=id_type, location=location, epochs=time2)

    if output == 'ephemeris':
        eph = ob.ephemerides(extra_precision=True)
        obstime = Time(eph['datetime_jd'], format='jd', scale='utc')
        pos = SkyCoord(eph['RA'], eph['DEC'], eph['delta'], frame='icrs', obstime=obstime)
    else:
        vec = ob.vectors(refplane='earth')
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
