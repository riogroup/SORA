import astropy.units as u
import numpy as np

__all__ = ['search_sbdb', 'apparent_magnitude']


def search_sbdb(name):
    """Searches JPL Small-Body DataBase to search for object information.

    As the name implies, it searches only for Small Bodies. Planets and
    satellites information are not retrieved by this function.

    Parameters
    ----------
    name : `str`
        The name of the object for the search. It can be the attributed `spkid`
        or `designation number`. The name is case insensitive.

    Returns
    -------
    sbdb : `dict`
        An ordered dictionary with the object information.

    Important
    ---------
    The query is not an autocomplete search, so ``name='Charikl'`` will not find
    Chariklo. If more than 1 object is found, the user is asked to select the
    correct one (e.g: ``name='neowise'``).
    """
    from astroquery.jplsbdb import SBDB
    from . import values

    print('Obtaining data for {} from SBDB'.format(name))
    sbdb = SBDB.query(name, full_precision=True, solution_epoch=True, validity=True, phys=True, discovery=True)
    if 'message' in sbdb:
        if sbdb['message'] == values.not_found_message:
            raise ValueError(values.not_found_message + " on SBDB")
        elif sbdb['message'] == values.many_objects_message:
            sbdb = select_body(sbdb)
    return sbdb


def select_body(sbdb):
    """Creates an object selection table.

    A table to select object is created when the SBDB search returns more than
    one object. This function is not supposed to be called by the user.

    Parameters
    ----------
    sbdb : `dict`
        An ordered dictionary object returned by SBDB search.

    Returns
    -------
    sbdb : `dict`
        An ordered dictionary with the data of the selected object.
    """
    from astropy.table import Table

    print('{} bodies were found.'.format(sbdb['count']))
    t = Table()
    t['num'] = np.arange(sbdb['count']) + 1
    t['Designation'] = sbdb['list']['pdes']
    t['Name'] = sbdb['list']['name']
    while True:
        t.pprint_all()
        print('0: Cancel')
        choice = int(input('Choose the corresponding number of the correct small body: '))
        if choice in np.arange(sbdb['count'] + 1):
            break
    if choice == 0:
        raise ValueError('It was not possible to define a Small Body')
    return search_sbdb(name=sbdb['list']['pdes'][choice - 1])


def apparent_magnitude(H, G, dist, sundist, phase=0.0):
    """Calculates the object's apparent magnitude.

    Parameters
    ----------
    H : `int`, `float`
        Absolute magnitude.

    G : `int`, `float`
        Slope parameter.

    dist : `int`, `float`
        Observer-object distance, in AU.

    sundist : `int`, `float`
        Sun-object distance, in AU.

    phase : `int`, `float`, default=0
        Phase angle (Sun-Target-Observer), in degrees.


    Returns
    -------
    ap_mag : `float`
        Apparent magnitude.
    """
    phi0 = np.exp(-3.33 * (np.tan(0.5 * phase * u.deg) ** 0.63))
    phi1 = np.exp(-1.87 * (np.tan(0.5 * phase * u.deg) ** 1.22))
    Ha = H - 2.5 * np.log10((1.0 - G) * phi0 + G * phi1)
    ap_mag = Ha + 5 * np.log10(sundist * dist)
    return ap_mag.value


# This block is a hard coded search for satellite information.
# For the moment it only covers some of what SBDB cannot offer.
# It is supposed to be replaced by a search on a trustful
# database on future versions.
grav = 'Grav et al. (2015). AJ 809(1):3'
rettig = 'Rettig et al. (2001). Icarus 154(2):313-320'
bauer = 'Bauer et al. (2006). Icarus 184(1):181-197'
veverka = 'Veverka et al. (1991). book Uranus 528-560'
satellites = {
    'io': {
        'spkid': 501,
        'albedo': [0.63, 0.02, 'Simonelli et al. (1984). Icarus 59:406-425'],
        'diameter': [3643.2, 1.0, 'horizons'],
        'GM': [5959.916, 0.012, 'JUP230']
    },
    'europa': {
        'spkid': 502,
        'albedo': [0.67, 0.03, 'Buratti et al. (1983). Icarus 55:93-110'],
        'diameter': [3121.6, 1.0, 'horizons'],
        'GM': [3202.739, 0.009, 'JUP230']
    },
    'ganymede': {
        'spkid': 503,
        'albedo': [0.43, 0.02, 'Morrison et al. (1977). book Planetary Satellites 363-378'],
        'diameter': [5262.4, 3.4, 'Anderson et al. (2000). BAAS 33(3):1101'],
        'GM': [9887.834, 0.017, 'JUP230']
    },
    'callisto': {
        'spkid': 504,
        'albedo': [0.17, 0.02, 'Morrison et al. (1977). book Planetary Satellites 363-378'],
        'diameter': [4820.6, 3.0, 'Morrison et al. (2000). Icarus 153:157-161'],
        'GM': [7179.289, 0.013, 'JUP230']
    },
    'amalthea': {
        'spkid': 505,
        'albedo': [0.090, 0.005, 'Simonelli et al. (2000). Icarus 147:353-365'],
        'diameter': [166.9, 4.8, 'horizons'],
        'GM': [0.138, 0.03, 'JUP230']
    },
    'thebe': {
        'spkid': 514,
        'albedo': [0.047, 0.003, 'Simonelli et al. (2000). Icarus 147:353-365'],
        'diameter': [98.6, 4.0, 'Thomas et al. (1998). Icarus 135:360-371'],
        'GM': [0.10, 0.0, 'horizons']
    },
    'himalia': {
        'spkid': 506,
        'albedo': [5.7, 0.1, grav],
        'H': [8.0, 0.01, rettig],
        'G': [0.10, 0.15, ''],
        'diameter': [139.6, 1.7, grav],
        'GM': [0.28, 0.04, 'Emelyanov (2005). A&A 438(3):L33-L36'],
        'rotation': [7.7819, 0.0005, 'Pilcher et al. (2012). Icarus 219(2):741-742']
    },
    'elara': {
        'spkid': 507,
        'albedo': [5.7, 0.8, grav],
        'H': [9.45, 0.01, rettig],
        'G': [0.10, 0.15, ''],
        'diameter': [79.9, 1.7, grav]
    },
    'pasiphae': {
        'spkid': 508,
        'albedo': [4.4, 0.6, grav],
        'H': [10.21, 0.01, rettig],
        'G': [0.10, 0.15, ''],
        'diameter': [57.8, 0.8, grav]
    },
    'sinope': {
        'spkid': 509,
        'albedo': [4.2, 0.6, grav],
        'H': [11.29, 0.01, rettig],
        'G': [0.10, 0.15, ''],
        'diameter': [35.0, 0.6, grav]
    },
    'lysithea': {
        'spkid': 510,
        'albedo': [3.6, 0.6, grav],
        'H': [11.09, 0.02, rettig],
        'G': [0.10, 0.15, ''],
        'diameter': [42.2, 0.7, grav]
    },
    'carme': {
        'spkid': 511,
        'albedo': [3.5, 0.6, grav],
        'H': [10.91, 0.01, rettig],
        'G': [0.10, 0.15, ''],
        'diameter': [46.7, 0.9, grav]
    },
    'ananke': {
        'spkid': 512,
        'albedo': [3.8, 0.6, grav],
        'H': [11.84, 0.03, rettig],
        'G': [0.10, 0.15, ''],
        'diameter': [29.1, 0.6, grav]
    },
    'leda': {
        'spkid': 513,
        'albedo': [3.4, 0.6, grav],
        'H': [12.63, 0.03, rettig],
        'G': [0.10, 0.15, ''],
        'diameter': [21.5, 1.7, grav]
    },
    'phoebe': {
        'spkid': 609,
        'albedo': [10.0, 0.5, grav],
        'H': [6.59, 0.02, bauer],
        'G': [0.02, 0.03, bauer],
        'diameter': [202.2, 4.5, grav],
        'density': [1.630, 0.045, 'Porco et al. (2005). Science 307(5713):1237-1242'],
        'GM': [0.5534, 0.0006, 'Jacobson et al. (2006). AJ 132(6):2520-2526'],
        'rotation': [9.27365, 0.00002, 'Gomes-Júnior et al. (2020). MNRAS 492(1):770-781'],
        'pole': [['356.6/77.9'], ['0/0'], 'Porco et al. (2005). Science 307(5713):1237-1242']
    },
    'albiorix': {
        'spkid': 626,
        'albedo': [6.2, 2.8, grav],
        'H': [11.35, 0.05, bauer],
        'G': [0.39, 0.06, bauer],
        'diameter': [28.6, 5.4, grav]
    },
    'siarnaq': {
        'spkid': 629,
        'albedo': [5.0, 1.7, grav],
        'H': [10.90, 0.05, bauer],
        'G': [0.45, 0.17, bauer],
        'diameter': [39.3, 5.9, grav]
    },
    'ariel': {
        'spkid': 701,
        'albedo': [0.39, 0.04, veverka],
        'H': [1.3, 0.2, veverka],
        'G': [0.10, 0.15, ''],
        'diameter': [1157.8, 1.2, 'Thomas et al. (1988). Icarus 73:427-441'],
        'GM': [86.4, 5.0, 'Jacobson et al. (2007). BAAS 39:453']
    },
    'umbriel': {
        'spkid': 702,
        'albedo': [0.21, 0.02, veverka],
        'H': [2.2, 0.2, veverka],
        'G': [0.10, 0.15, ''],
        'diameter': [1169.4, 5.6, 'Thomas et al. (1988). Icarus 73:427-441'],
        'GM': [81.5, 5.0, 'Jacobson et al. (2007). BAAS 39:453']
    },
    'titania': {
        'spkid': 703,
        'albedo': [0.27, 0.03, veverka],
        'H': [1.2, 0.2, veverka],
        'G': [0.10, 0.15, ''],
        'diameter': [1577.8, 3.6, 'Thomas et al. (1988). Icarus 73:427-441'],
        'GM': [228.2, 5.0, 'Jacobson et al. (2007). BAAS 39:453']
    },
    'oberon': {
        'spkid': 704,
        'albedo': [0.23, 0.03, veverka],
        'H': [1.4, 0.2, veverka],
        'G': [0.10, 0.15, ''],
        'diameter': [1522.8, 5.2, 'Thomas et al. (1988). Icarus 73:427-441'],
        'GM': [192.4, 7.0, 'Jacobson et al. (2007). BAAS 39:453']
    },
    'miranda': {
        'spkid': 705,
        'albedo': [0.32, 0.03, veverka],
        'H': [3.4, 0.0, veverka],
        'G': [0.10, 0.15, ''],
        'diameter': [471.6, 1.4, 'Thomas et al. (1988). Icarus 73:427-441'],
        'GM': [4.4, 0.4, 'Jacobson et al. (2007). BAAS 39:453']
    },
    'caliban': {
        'spkid': 716,
        'albedo': [0.04, 0.0, 'Sheppard et al. (2005). AJ 129(1):518-525'],
        'H': [9.16, 0.04, 'Grav et al. (2004). AJ 613:L77–L80'],
        'G': [0.10, 0.15, ''],
        'diameter': [72, 0.0, 'Sheppard et al. (2005). AJ 129(1):518-525'],
        'rotation': [6.9190, 0.0082, 'Farkas-Takács et al. (2017) 154(3):119']
    },
    'sycorax': {
        'spkid': 717,
        'albedo': [0.04, 0.0, 'Sheppard et al. (2005). AJ 129(1):518-525'],
        'H': [7.50, 0.04, 'Grav et al. (2004). AJ 613:L77–L80'],
        'G': [0.10, 0.15, ''],
        'diameter': [150, 0.0, 'Sheppard et al. (2005). AJ 129(1):518-525'],
        'rotation': [9.948, 0.019, 'Farkas-Takács et al. (2017) 154(3):119']
    },
    'nereid': {
        'spkid': 802,
        'albedo': [0.155, 0.0, 'Thomas et al. (1991). JGR 96(19):253'],
        'H': [4.44, 0.01, 'Grav et al. (2004). AJ 613:L77–L80'],
        'G': [0.10, 0.15, ''],
        'diameter': [340, 50, 'Thomas et al. (1991). JGR 96(19):253']
    }
}


def search_satdb(name):
    sat = satellites.get(name.lower())
    if sat is None:
        raise ValueError('specified object was not found')
    return sat
# end of hard coded block.
