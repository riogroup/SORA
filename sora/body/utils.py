from . import values
from astroquery.jplsbdb import SBDB
import numpy as np
import astropy.units as u
from astropy.table import Table


__all__ = ['search_sbdb', 'apparent_magnitude']


def search_sbdb(name):
    """Searchs JPL Small-Body DataBase to search for object information.
    As the name implies, it looks only for Small Bodies.
    Planets and satellites information are not retrieved by this function.

    Parameters:
        name (str): The name of the object for the search.
            It can be the used spkid or designation number.
            The name is case insensitive.

    Return:
        sbdb (dict): An OrderedDict with the object information

    Important note:
        The query is not an autocomplete search, so name='Charikl'
        will not find Chariklo. If more than 1 object is found, the user
        is asked to select the correct one (e.g: name='neowise')
    """
    print('Obtaining data for {} from SBDB'.format(name))
    sbdb = SBDB.query(name, full_precision=True, solution_epoch=True, validity=True, phys=True, discovery=True)
    if 'message' in sbdb:
        if sbdb['message'] == values.not_found_message:
            raise ValueError(values.not_found_message + " on SBDB")
        elif sbdb['message'] == values.many_objects_message:
            sbdb = select_body(sbdb)
    return sbdb


def select_body(sbdb):
    """Creates table to select object when SBDB search returns more than 1 object.
    This function is not supposed to be called by the user.

    Parameters:
        sbdb (dict): An OrderedDict object returned by SBDB search.

    Return:
        sbdb (dict): An OrderedDict with the data of the selected object.
    """
    print('{} bodies were found.'.format(sbdb['count']))
    t = Table()
    t['num'] = np.arange(sbdb['count'])+1
    t['Designation'] = sbdb['list']['pdes']
    t['Name'] = sbdb['list']['name']
    while True:
        t.pprint_all()
        print('0: Cancel')
        choice = int(input('Choose the corresponding number of the correct small body: '))
        if choice in np.arange(sbdb['count']+1):
            break
    if choice == 0:
        raise ValueError('It was not possible to define a Small Body')
    return search_sbdb(name=sbdb['list']['pdes'][choice-1])


def apparent_magnitude(H, G, dist, sundist, phase=0.0):
    """ Calculates the Object Apparent Magnitude

    Parameters:
        H (int, float): Absolute Magnitude
        G (int, float): Slope parameter
        dist (int, float): Observer-Object distance, in AU
        sundist (int, float): Sun-Object distance, in AU
        phase (int, float): Phase Angle: Sun-Target-Observer, in deg

    Returns:
        ap_mag (float): Apparent Magnitude
    """
    phi0 = np.exp(-3.33*(np.tan(0.5*phase*u.deg)**0.63))
    phi1 = np.exp(-1.87*(np.tan(0.5*phase*u.deg)**1.22))
    Ha = H - 2.5*np.log10((1.0-G)*phi0 + G*phi1)
    ap_mag = Ha + 5*np.log10(sundist*dist)
    return ap_mag.value
