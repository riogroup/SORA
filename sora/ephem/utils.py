import warnings

warnings.simplefilter('always', UserWarning)

__all__ = ['getBSPfromJPL', 'ephem_kernel']


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


def ephem_kernel(time, target, observer, kernels):
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

    if type(kernels) == str:
        kernels = [kernels]
    for kern in kernels:
        spice.furnsh(kern)
    time = Time(time)
    t0 = Time('J2000', scale='tdb')
    if time.isscalar:
        time = Time([time])
    dt = (time - t0)
    delt = 0 * u.s
    # calculates vector Observer -> Solar System Barycenter
    position1 = np.array(spice.spkpos('0', dt.sec, 'J2000', 'NONE', observer)[0])
    while True:
        # calculates new time
        tempo = dt - delt
        # calculates vector Solar System Barycenter -> Object
        position2 = spice.spkpos(target, tempo.sec, 'J2000', 'NONE', '0')[0]
        position = (position1 + position2).T
        # calculates linear distance Earth Topocenter -> Object
        dist = np.linalg.norm(position, axis=0) * u.km
        # calculates new light time
        delt = (dist / const.c).decompose()
        # if difference between new and previous light time is smaller than 0.001 sec, than continue.
        if all(np.absolute(((dt - tempo) - delt).sec) < 0.001):
            break
    coord = SkyCoord(position[0], position[1], position[2], frame='icrs', unit=u.km,
                     representation_type='cartesian', obstime=time)
    spice.kclear()
    coord_rd = SkyCoord(ra=coord.spherical.lon, dec=coord.spherical.lat,
                        distance=coord.spherical.distance, obstime=time)
    if len(coord) == 1:
        return coord_rd[0]
    return coord_rd
