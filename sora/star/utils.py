from astropy.coordinates import SkyCoord
from astropy.table import Table
import astropy.units as u
import astropy.constants as const
from astroquery.vizier import Vizier
import numpy as np
from sora.config import input_tests


__all__ = ['van_belle', 'kervella']


def search_star(**kwargs):
    """ Searches position on VizieR and returns a catalogue.

    Parameters:
        coord (str, SkyCoord): Coordinate to perform the search.
        code (str): Gaia Source_id of the star
        columns (list): List of strings with the name of the columns to retrieve.
        radius (int, float, unit.quantity): Radius to search around coordinate
        catalog (str): VizieR catalogue to search.

    Returns:
        catalogue(astropy.Table): An astropy Table with the catalogue informations.
    """
    input_tests.check_kwargs(kwargs, allowed_kwargs=['catalog', 'code', 'columns', 'coord', 'log', 'radius'])
    row_limit = 100
    if 'log' in kwargs and kwargs['log']:
        print('\nDownloading star parameters from {}'.format(kwargs['catalog']))
    vquery = Vizier(columns=kwargs['columns'], row_limit=row_limit, timeout=600)
    if 'code' in kwargs:
        catalogue = vquery.query_constraints(catalog=kwargs['catalog'], Source=kwargs['code'], cache=False)
    elif 'coord' in kwargs:
        catalogue = vquery.query_region(kwargs['coord'], radius=kwargs['radius'], catalog=kwargs['catalog'], cache=False)
    else:
        raise ValueError('At least a code or coord should be given as input')
    return catalogue


def van_belle(magB=None, magV=None, magK=None):
    """ Determines the diameter of a star in mas using equations from van Belle (1999)
        -- Publi. Astron. Soc. Pacific 111, 1515-1523:

    Parameters:
        magB: The magnitude B of the star
        magV: The magnitude V of the star
        magK: The magnitude K of the star
        If any of those values is 'None', 'nan' or higher than 49, it is not considered.
    """
    if magB is None or np.isnan(magB) or magB > 49:
        magB = np.nan
    if magV is None or np.isnan(magV) or magV > 49:
        magV = np.nan
    if magK is None or np.isnan(magK) or magK > 49:
        magK = np.nan

    def calc_diameter(a1, a2, mag):
        return 10**(a1 + a2*(mag - magK) - 0.2*mag)

    params = {'sg': {'B': [0.648, 0.220], 'V': [0.669, 0.223]},
              'ms': {'B': [0.500, 0.290], 'V': [0.500, 0.264]},
              'vs': {'B': [0.840, 0.211], 'V': [0.789, 0.218]}}

    mag = np.array([magB, magV])
    diameter = {}
    for st in ['sg', 'ms', 'vs']:
        diameter_s = {}
        for i, m in enumerate(['B', 'V']):
            diam = calc_diameter(*params[st][m], mag[i])
            if not np.isnan(diam):
                diameter_s[m] = calc_diameter(*params[st][m], mag[i])*u.mas
        if diameter_s:
            diameter[st] = diameter_s
    return diameter


def kervella(magB=None, magV=None, magK=None):
    """ Determines the diameter of a star in mas using equations from Kervella et. al (2004)
        -- A&A Vol. 426, No.  1:

    Parameters:
        magB: The magnitude B of the star
        magV: The magnitude V of the star
        magK: The magnitudes K of the star
        If any of those values is 'None', 'nan' or higher than 49, it is not considered.
    """
    if magB is None or np.isnan(magB) or magB > 49:
        magB = np.nan
    if magV is None or np.isnan(magV) or magV > 49:
        magV = np.nan
    if magK is None or np.isnan(magK) or magK > 49:
        magK = np.nan
    const1 = np.array([0.0755, 0.0535])
    const2 = np.array([0.5170, 0.5159])
    mag = np.array([magV, magB])
    vals = 10**(const1*(mag-magK)+const2-0.2*magK)
    diam = {}
    if not np.isnan(vals[0]):
        diam['V'] = vals[0]*u.mas
    if not np.isnan(vals[1]):
        diam['B'] = vals[1]*u.mas
    return diam


def spatial_motion(ra, dec, pmra, pmdec, parallax=0, rad_vel=0,  dt=0, cov_matrix=None):
    """ Applies spatial motion to star coordinate

    Parameters:
        ra (int, float): Right Ascension of the star at t=0 epoch, in deg.
        dec (int, float): Declination of the star at t=0 epoch, in deg.
        pmra (int, float): Proper Motion in RA of the star at t=0 epoch, in mas/year.
        pmdec (int, float): Proper Motion in DEC of the star at t=0 epoch, in mas/year.
        parallax (int, float): Parallax of the star at t=0 epoch, in mas.
        rad_vel (int, float): Radial Velocity of the star at t=0 epoch, in km/s.
        dt (int, float): Variation of time from catalogue epoch, in days.
        cov_matrix (2D-array): 6x6 covariance matrix.
    """
    A = (1*u.AU).to(u.km).value  # Astronomical units in km
    c = const.c.to(u.km/u.year).value  # light velocity

    par = True
    # Eliminate negative or zero parallaxes
    if parallax is None or parallax <= 0:
        par = False
        parallax = 1e-4

    if cov_matrix is not None and cov_matrix.shape != (6, 6):
        raise ValueError('Covariance matrix must be a 6x6 matrix')

    ra0 = u.Quantity(ra, unit=u.deg).to(u.rad).value
    dec0 = u.Quantity(dec, unit=u.deg).to(u.rad).value
    parallax0 = u.Quantity(parallax, unit=u.mas).to(u.rad).value
    pmra0 = u.Quantity(pmra, unit=u.mas/u.year).to(u.rad/u.year).value
    pmdec0 = u.Quantity(pmdec, unit=u.mas/u.year).to(u.rad/u.year).value
    rad_vel0 = u.Quantity(rad_vel, unit=u.km/u.s).to(u.AU/u.year).value
    dt = u.Quantity(dt, unit=u.day).to(u.year).value

    # normal triad relative to the celestial sphere
    # p0 points to growing RA, q0 to growing DEC and r0 to growing distance.
    p0 = np.array([-np.sin(ra0), np.cos(ra0), 0.0]).T
    q0 = np.array([-np.sin(dec0)*np.cos(ra0), -np.sin(dec0)*np.sin(ra0), np.cos(dec0)])
    r0 = np.array([np.cos(dec0)*np.cos(ra0), np.cos(dec0)*np.sin(ra0), np.sin(dec0)])

    b0 = A/parallax0
    tau_0 = b0/c
    tau_A = A/c

    vec_b0 = b0*r0  # distance vector
    vec_u0 = vec_b0/np.linalg.norm(vec_b0)
    vec_mi0 = np.array(p0*pmra0 + q0*pmdec0)  # proper motion vector

    mi_r0 = rad_vel0/b0
    mi0 = np.sqrt(pmra0**2+pmdec0**2)  # total proper motion

    v0 = b0*(r0*mi_r0+vec_mi0)  # apparent space velocity
    v_r0 = np.linalg.norm(v0)

    # Scaling factors of time, distance and velocity due to light time
    f_T = ((dt + 2*tau_0)/(tau_0+(1-v_r0/c)*dt + np.sqrt(np.linalg.norm((vec_b0+v0*dt))**2
           + (2*dt/(c**2*tau_0))*np.linalg.norm(np.cross(v0, vec_b0))**2)/c))
    f_D = np.sqrt(1+2*mi_r0*dt*f_T + (mi0**2 + mi_r0**2)*(dt*f_T)**2)
    f_V = (1 + (tau_A/parallax0)*(mi_r0*(f_D - 1) + f_D*(mi0**2 + mi_r0**2)*dt*f_T))

    vec_u = (r0*(1 + mi_r0*dt*f_T) + vec_mi0*dt*f_T)*f_D
    vec_mi = (vec_mi0*(1 + mi_r0*dt*f_T) - vec_u0*mi0**2*f_T)*f_D**3*f_V
    mi_r = (mi_r0 + (mi0**2 + mi_r0**2)*dt*f_T)*f_D**2*f_V

    dec = np.arcsin(vec_u[2])  # new dec
    ra = np.arctan2(vec_u[1]/np.cos(dec), vec_u[0]/np.cos(dec))  # new ra

    parallax = parallax0*f_D  # new parallax
    new_dist = A/parallax  # new distance

    if par:
        coord = SkyCoord(ra*u.rad, dec*u.rad, new_dist*u.km)
    else:
        coord = SkyCoord(ra*u.rad, dec*u.rad)

    if cov_matrix is None:
        return coord

    p = np.array([-np.sin(ra), np.cos(ra), 0.0])
    q = np.array([-np.sin(dec)*np.cos(ra), -np.sin(dec)*np.sin(ra), np.cos(dec)])

    Z = np.sqrt(1 + (dt + 2*tau_A/parallax0)*mi0**2*dt + (2 + mi_r0*dt)*mi_r0*dt)
    Y = parallax0*dt + tau_A*(1 + Z - mi_r0*dt)
    X = parallax0*dt + 2*tau_A

    # partial derivatives of the logarithm of the velocity factor
    # chi_parallax = (1/parallax0)*(1 - f_V)  # not used
    chi_pm = (tau_A/parallax0)*dt*f_T*f_D*(mi_r*dt*f_T - 2*f_V)
    chi_r = (tau_A/parallax0)*(f_V + f_D*(f_V + (1 + mi_r0*dt*f_T)*(mi_r*dt*f_T - 2*f_V)))
    # chi_T = -(tau_A/parallax0)*f_D**3*mi0**2*dt*f_T  # not used

    # partial derivatives of the logarithmic of the time factor
    psi_parallax = dt/X - (dt/Y)*(1 - (mi0**2*tau_A**2)/(Z*parallax0**2))
    psi_pm = -((dt*tau_A)/(Y*Z))*(dt + 2*tau_A/parallax0)
    psi_r = -(dt*tau_A/Y)*((1 + mi_r0*dt)/Z - 1)

    ni = vec_mi*(1 - dt*f_T*(3*mi_r/f_V + (tau_A/parallax0)*mi0**2*f_D**3*f_V)) - vec_mi0*f_D**3*f_V
    eta = mi_r*(1 - dt*f_T*(2*mi_r/f_V + (tau_A/parallax0)*mi0**2*f_D**3*f_V)) - mi_r0*f_D**2*f_V

    p_l = p/np.linalg.norm(p)
    q_l = q/np.linalg.norm(q)

    pmra = np.dot(p_l, vec_mi)  # new pmra
    pmdec = np.dot(q_l, vec_mi)  # new pmdec

    # Jacobian matrix
    J = np.zeros((6, 6))

    # d(alpha)/d(valores inicias)
    J[0, 0] = np.dot(p_l, p0)*(1 + mi_r0*dt*f_T)*f_D - np.dot(p_l, r0)*pmra0*dt*f_T*f_D
    J[0, 1] = np.dot(p_l, q0)*(1 + mi_r0*dt*f_T)*f_D - np.dot(p_l, r0)*pmdec0*dt*f_T*f_D
    J[0, 2] = - np.dot(p_l, r0)*f_D*psi_parallax
    J[0, 3] = np.dot(p_l, p0)*dt*f_T*f_D - np.dot(p_l, r0)*pmra0*f_D*psi_pm
    J[0, 4] = np.dot(p_l, q0)*dt*f_T*f_D - np.dot(p_l, r0)*pmdec0*f_D*psi_pm
    J[0, 5] = - pmra*(dt*f_T)**2/f_V - np.dot(p_l, r0)*f_D*psi_r

    # d(delta)/d(valores inicias)
    J[1, 0] = np.dot(q_l, p0)*(1 + mi_r0*dt*f_T)*f_D - np.dot(q_l, r0)*pmra0*dt*f_T*f_D
    J[1, 1] = np.dot(q_l, q0)*(1 + mi_r0*dt*f_T)*f_D - np.dot(q_l, r0)*pmdec0*dt*f_T*f_D
    J[1, 2] = - np.dot(q_l, r0)*f_D*psi_parallax
    J[1, 3] = np.dot(q_l, p0)*dt*f_T*f_D - np.dot(q_l, r0)*pmra0*f_D*psi_pm
    J[1, 4] = np.dot(q_l, q0)*dt*f_T*f_D - np.dot(q_l, r0)*pmdec0*f_D*psi_pm
    J[1, 5] = - pmdec*(dt*f_T)**2/f_V - np.dot(q_l, r0)*f_D*psi_r

    # d(parallax)/d(valores inicias)
    J[2, 0] = 0
    J[2, 1] = 0
    J[2, 2] = f_D - parallax*(mi_r*dt*f_T/f_V)*psi_parallax
    J[2, 3] = - parallax*pmra0*((dt*f_T)**2*f_D**2 + (mi_r*dt*f_T/f_V)*psi_pm)
    J[2, 4] = - parallax*pmdec0*((dt*f_T)**2*f_D**2 + (mi_r*dt*f_T/f_V)*psi_pm)
    J[2, 5] = - parallax*((1 + mi_r0*dt*f_T)*dt*f_T*f_D**2 + (mi_r*dt*f_T/f_V)*psi_r)

    # d(pmra)/d(valores inicias)
    J[3, 0] = - np.dot(p_l, p0)*mi0**2*dt*f_T*f_D**3*f_V - np.dot(p_l, r0)*pmra0*(1+mi_r0*dt*f_T)*f_D**3*f_V
    J[3, 1] = - np.dot(p_l, q0)*mi0**2*dt*f_T*f_D**3*f_V - np.dot(p_l, r0)*pmdec0*(1+mi_r0*dt*f_T)*f_D**3*f_V
    J[3, 2] = np.dot(p_l, ni)*psi_parallax
    J[3, 3] = np.dot(p_l, p0)*(1 + mi_r0*dt*f_T)*f_D**3*f_V - 2*np.dot(p_l, r0)*pmra0*dt*f_T*f_D**3*f_V \
        - 3*pmra*pmra0*(dt*f_T)**2*f_D**2*f_V + pmra*pmra0*chi_pm + np.dot(p_l, ni)*pmra0*psi_pm
    J[3, 4] = np.dot(p_l, q0)*(1 + mi_r0*dt*f_T)*f_D**3*f_V - 2*np.dot(p_l, r0)*pmdec0*dt*f_T*f_D**3*f_V \
        - 3*pmra*pmdec0*(dt*f_T)**2*f_D**2*f_V + pmra*pmdec0*chi_pm + np.dot(p_l, ni)*pmdec0*psi_pm
    J[3, 5] = np.dot(p_l, (vec_mi0*f_D - 3*vec_mi*(1 + mi_r0*dt*f_T)))*dt*f_T*f_D**2*f_V

    # d(pmdec)/d(valores inicias)
    J[4, 0] = - np.dot(q_l, p0)*mi0**2*dt*f_T*f_D**3*f_V - np.dot(q_l, r0)*pmra0*(1+mi_r0*dt*f_T)*f_D**3*f_V
    J[4, 1] = - np.dot(q_l, q0)*mi0**2*dt*f_T*f_D**3*f_V - np.dot(q_l, r0)*pmdec0*(1+mi_r0*dt*f_T)*f_D**3*f_V
    J[4, 2] = np.dot(q_l, ni)*psi_parallax
    J[4, 3] = np.dot(q_l, p0)*(1 + mi_r0*dt*f_T)*f_D**3*f_V - 2*np.dot(q_l, r0)*pmra0*dt*f_T*f_D**3*f_V \
        - 3*pmdec*pmra0*(dt*f_T)**2*f_D**2*f_V + pmdec*pmra0*chi_pm + np.dot(q_l, ni)*pmra0*psi_pm
    J[4, 4] = np.dot(q_l, q0)*(1 + mi_r0*dt*f_T)*f_D**3*f_V - 2*np.dot(q_l, r0)*pmdec0*dt*f_T*f_D**3*f_V \
        - 3*pmdec*pmdec0*(dt*f_T)**2*f_D**2*f_V + pmdec*pmdec0*chi_pm + np.dot(q_l, ni)*pmdec0*psi_pm
    J[4, 5] = np.dot(q_l, (vec_mi0*f_D - 3*vec_mi*(1 + mi_r0*dt*f_T)))*dt*f_T*f_D**2*f_V

    # d(rad_vel)/d(valores inicias)
    J[5, 0] = 0
    J[5, 1] = 0
    J[5, 2] = eta*psi_parallax
    J[5, 3] = 2*pmra0*(1 + mi_r0*dt*f_T)*dt*f_T*f_D**4*f_V + pmra0*mi_r*chi_pm + pmra0*eta*psi_pm
    J[5, 4] = 2*pmdec0*(1 + mi_r0*dt*f_T)*dt*f_T*f_D**4*f_V + pmdec0*mi_r*chi_pm + pmdec0*eta*psi_pm
    J[5, 5] = ((1 + mi_r0*dt*f_T)**2 - mi0**2*(dt*f_T)**2)*f_D**4*f_V + mi_r*chi_r + eta*psi_r

    # Propagated covariance matrix
    C = np.matmul(J, np.matmul(cov_matrix, J.T))

    err = np.array([np.sqrt(C[i, i]) for i in np.arange(3)])
    if not par:
        err = err[:2]

    return coord, err


def choice_star(catalogue, coord, columns, source):
    tstars = SkyCoord(catalogue[columns[0]], catalogue[columns[1]])
    sep = tstars.separation(coord)
    k = sep.argsort()
    while True:
        t = Table()
        t['num'] = np.arange(len(tstars))+1
        t['dist(")'] = sep[k].arcsec
        t['dist(")'].format = '6.3f'
        for c in columns[2:]:
            t[c] = catalogue[c][k].quantity.value
            t[c].format = '6.3f'
        t['RA___ICRS___DEC'] = tstars[k].to_string('hmsdms', precision=4)
        t.pprint_all()
        print('  0: Cancel')
        choice = int(input('Choose the corresponding number of the correct star: '))
        if choice in np.arange(len(k)+1):
            break
        print('{} is not a valid choice. Please select the correct star'.format(choice))
    if choice == 0:
        if source == 'gaia':
            raise ValueError('It was not possible to define a star')
        elif source == 'nomad':
            print('No magnitudes were obtained from NOMAD')
            return
        elif source == 'bjones':
            print('It was not possible to define a star')
            return
    return catalogue[[k[choice-1]]]
