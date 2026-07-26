import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord

from sora.config import input_tests
from sora.config.decorators import deprecated_alias

__all__ = ['van_belle', 'kervella']


@deprecated_alias(log='verbose')  # remove this line for v1.0
def search_star(**kwargs):
    """Searches a star position on VizieR and returns a catalogue result.

    Parameters
    ----------
    coord : `str`, `astropy.coordinates.SkyCoord`
        Coordinate to perform the search.

    code : `str`
        Gaia source identifier of the star.

    columns : `list`
        List of strings with the name of the columns to retrieve.

    radius : `int`, `float`, `astropy.units.Quantity`
        Radius to search around coordinates.

    catalog : `str`
        VizieR catalogue to search.

    verbose : `bool`
        If True, prints the catalogue being queried.

    Returns
    -------
    catalogue : `astroquery.utils.commons.TableList`
        Query result with catalogue information.
    """
    from astroquery.vizier import Vizier

    input_tests.check_kwargs(kwargs, allowed_kwargs=['catalog', 'code', 'columns', 'coord', 'verbose', 'radius'])
    row_limit = 100
    if 'verbose' in kwargs and kwargs['verbose']:
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
    """Determines the diameter of a star in mas using equations from van Belle (1999).

    See: Publications of the Astronomical Society of the Pacific, 111, 1515-1523.

    Parameters
    ----------
    magB : `float`, default=None
        The magnitude B of the star.

    magV : `float`, default=None
        The magnitude V of the star.

    magK : `float`, default=None
        The magnitude K of the star.


    Returns
    -------
    diameter : `dict`
        Angular diameters by stellar type and band.

    Notes
    -----
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
    """Determines the diameter of a star in mas using equations from Kervella et al. (2004).

    See: Astronomy & Astrophysics, 426, 297-307.

    Parameters
    ----------
    magB : `float`, default=None
        The magnitude B of the star.

    magV : `float`, default=None
        The magnitude V of the star.

    magK : `float`, default=None
        The magnitude K of the star.

    Returns
    -------
    diameter : `dict`
        Angular diameters by band.

    Notes
    -----
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
    diameter = {}
    if not np.isnan(vals[0]):
        diameter['V'] = vals[0]*u.mas
    if not np.isnan(vals[1]):
        diameter['B'] = vals[1]*u.mas
    return diameter


def spatial_motion(ra, dec, pmra, pmdec, parallax=0, rad_vel=0,  dt=0, cov_matrix=None):
    """Applies spatial motion to a star coordinate.

    This function supports either one star propagated to several ``dt`` values,
    or several stars propagated to one shared ``dt`` value. Astrometric
    parameter arrays are paired item by item and must all have the same shape;
    no broadcast is performed between star parameters. If ``cov_matrix`` is
    provided, only the one-star, one-or-many-``dt`` case is supported.

    Parameters
    ----------
    ra : `int`, `float`, `astropy.units.Quantity`
        Right Ascension of the star at t=0 epoch, in deg.

    dec : `int`, `float`, `astropy.units.Quantity`
        Declination of the star at t=0 epoch, in deg.

    pmra : `int`, `float`, `astropy.units.Quantity`
        Proper motion in RA*cos(DEC) of the star at t=0 epoch, in mas/year.

    pmdec : `int`, `float`, `astropy.units.Quantity`
        Proper motion in DEC of the star at t=0 epoch, in mas/year.

    parallax : `int`, `float`, `astropy.units.Quantity`, default=0
        Parallax of the star at t=0 epoch, in mas.

    rad_vel : `int`, `float`, `astropy.units.Quantity`, default=0
        Radial velocity of the star at t=0 epoch, in km/s.

    dt : `int`, `float`, `astropy.units.Quantity`, default=0
        Variation of time from catalogue epoch, in days.

    cov_matrix : `2D-array`, optional
        6x6 covariance matrix. It can only be used with one star. When ``dt``
        has multiple values, one propagated error is returned for each value.

    Returns
    -------
    coord : `astropy.coordinates.SkyCoord`
        Star coordinate propagated by spatial motion.
    errors : `numpy.ndarray`
        Returned only when ``cov_matrix`` is provided. Propagated coordinate
        uncertainties.
    """
    import astropy.constants as const

    A = (1*u.AU).to(u.km).value  # Astronomical units in km
    c = const.c.to(u.km/u.year).value  # light velocity

    if parallax is None:
        parallax = 0

    ra0 = np.asarray(u.Quantity(ra, unit=u.deg).to(u.rad).value, dtype=float)
    dec0 = np.asarray(u.Quantity(dec, unit=u.deg).to(u.rad).value, dtype=float)
    parallax_mas = np.asarray(u.Quantity(parallax, unit=u.mas).to(u.mas).value, dtype=float)
    pmra0 = np.asarray(u.Quantity(pmra, unit=u.mas/u.year).to(u.rad/u.year).value, dtype=float)
    pmdec0 = np.asarray(u.Quantity(pmdec, unit=u.mas/u.year).to(u.rad/u.year).value, dtype=float)
    rad_vel0 = np.asarray(u.Quantity(rad_vel, unit=u.km/u.s).to(u.AU/u.year).value, dtype=float)
    dt = np.asarray(u.Quantity(dt, unit=u.day).to(u.year).value, dtype=float)

    star_shapes = {
        'ra': np.shape(ra0),
        'dec': np.shape(dec0),
        'pmra': np.shape(pmra0),
        'pmdec': np.shape(pmdec0),
        'parallax': np.shape(parallax_mas),
        'rad_vel': np.shape(rad_vel0),
    }
    non_scalar_shapes = [shape for shape in star_shapes.values() if shape != ()]
    star_shape = non_scalar_shapes[0] if non_scalar_shapes else ()
    if any(shape != star_shape for shape in star_shapes.values()):
        raise ValueError(
            'Star astrometric parameters (ra, dec, pmra, pmdec, parallax, rad_vel) '
            'must all be scalars or have exactly the same shape.'
        )

    star_size = int(np.prod(star_shape)) if star_shape else 1
    dt_size = int(np.prod(dt.shape)) if dt.shape else 1
    if star_size > 1 and dt_size > 1:
        raise ValueError('spatial_motion accepts either several stars and one dt, or one star and several dt values.')

    if cov_matrix is not None:
        cov_matrix = np.asarray(cov_matrix, dtype=float)
        if cov_matrix.shape != (6, 6):
            raise ValueError('Covariance matrix must be a 6x6 matrix')
        if star_size > 1:
            raise ValueError('Covariance matrix propagation is only supported for one star.')

    # Eliminate negative, zero, or invalid parallaxes for the propagation itself.
    par = np.isfinite(parallax_mas) & (parallax_mas > 0)
    parallax0 = u.Quantity(np.where(par, parallax_mas, 1e-4), unit=u.mas).to(u.rad).value

    # normal triad relative to the celestial sphere
    # p0 points to growing RA, q0 to growing DEC and r0 to growing distance.
    p0 = np.stack([-np.sin(ra0), np.cos(ra0), np.zeros_like(ra0)], axis=-1)
    q0 = np.stack([-np.sin(dec0)*np.cos(ra0), -np.sin(dec0)*np.sin(ra0), np.cos(dec0)], axis=-1)
    r0 = np.stack([np.cos(dec0)*np.cos(ra0), np.cos(dec0)*np.sin(ra0), np.sin(dec0)], axis=-1)

    b0 = A/parallax0
    tau_0 = b0/c
    tau_A = A/c

    vec_b0 = b0[..., None]*r0  # distance vector
    vec_u0 = vec_b0/np.linalg.norm(vec_b0, axis=-1, keepdims=True)
    vec_mi0 = p0*pmra0[..., None] + q0*pmdec0[..., None]  # proper motion vector

    mi_r0 = rad_vel0/b0
    mi0 = np.sqrt(pmra0**2+pmdec0**2)  # total proper motion

    v0 = b0[..., None]*(r0*mi_r0[..., None]+vec_mi0)  # apparent space velocity
    v_r0 = np.linalg.norm(v0, axis=-1)

    # Scaling factors of time, distance and velocity due to light time
    v0_dt = vec_b0 + v0*dt[..., None]
    cross_norm = np.linalg.norm(np.cross(v0, vec_b0), axis=-1)
    f_T = ((dt + 2*tau_0)/(tau_0+(1-v_r0/c)*dt + np.sqrt(np.linalg.norm(v0_dt, axis=-1)**2
           + (2*dt/(c**2*tau_0))*cross_norm**2)/c))
    f_D = np.sqrt(1+2*mi_r0*dt*f_T + (mi0**2 + mi_r0**2)*(dt*f_T)**2)
    f_V = (1 + (tau_A/parallax0)*(mi_r0*(f_D - 1) + f_D*(mi0**2 + mi_r0**2)*dt*f_T))

    vec_u = (r0*(1 + mi_r0*dt*f_T)[..., None] + vec_mi0*(dt*f_T)[..., None])*f_D[..., None]
    vec_mi = (vec_mi0*(1 + mi_r0*dt*f_T)[..., None] - vec_u0*(mi0**2*f_T)[..., None])*(f_D**3*f_V)[..., None]
    mi_r = (mi_r0 + (mi0**2 + mi_r0**2)*dt*f_T)*f_D**2*f_V

    dec = np.arcsin(vec_u[..., 2])  # new dec
    ra = np.arctan2(vec_u[..., 1]/np.cos(dec), vec_u[..., 0]/np.cos(dec))  # new ra

    parallax = parallax0*f_D  # new parallax
    new_dist = A/parallax  # new distance

    par_out = np.broadcast_to(par, np.shape(ra))
    if np.all(par_out):
        coord = SkyCoord(ra*u.rad, dec*u.rad, new_dist*u.km)
    elif not np.any(par_out):
        coord = SkyCoord(ra*u.rad, dec*u.rad)
    else:
        coord = SkyCoord(ra*u.rad, dec*u.rad, np.where(par_out, new_dist, np.nan)*u.km)

    if cov_matrix is None:
        return coord

    def dot_last(left, right):
        return np.sum(left*right, axis=-1)

    p = np.stack([-np.sin(ra), np.cos(ra), np.zeros_like(ra)], axis=-1)
    q = np.stack([-np.sin(dec)*np.cos(ra), -np.sin(dec)*np.sin(ra), np.cos(dec)], axis=-1)

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

    ni = vec_mi*(1 - dt*f_T*(3*mi_r/f_V + (tau_A/parallax0)*mi0**2*f_D**3*f_V))[..., None] - \
        vec_mi0*(f_D**3*f_V)[..., None]
    eta = mi_r*(1 - dt*f_T*(2*mi_r/f_V + (tau_A/parallax0)*mi0**2*f_D**3*f_V)) - mi_r0*f_D**2*f_V

    p_l = p/np.linalg.norm(p, axis=-1, keepdims=True)
    q_l = q/np.linalg.norm(q, axis=-1, keepdims=True)

    pmra = dot_last(p_l, vec_mi)  # new pmra
    pmdec = dot_last(q_l, vec_mi)  # new pmdec

    # Jacobian matrix
    J = np.zeros(np.shape(ra) + (6, 6))

    # d(alpha)/d(valores inicias)
    J[..., 0, 0] = dot_last(p_l, p0)*(1 + mi_r0*dt*f_T)*f_D - dot_last(p_l, r0)*pmra0*dt*f_T*f_D
    J[..., 0, 1] = dot_last(p_l, q0)*(1 + mi_r0*dt*f_T)*f_D - dot_last(p_l, r0)*pmdec0*dt*f_T*f_D
    J[..., 0, 2] = - dot_last(p_l, r0)*f_D*psi_parallax
    J[..., 0, 3] = dot_last(p_l, p0)*dt*f_T*f_D - dot_last(p_l, r0)*pmra0*f_D*psi_pm
    J[..., 0, 4] = dot_last(p_l, q0)*dt*f_T*f_D - dot_last(p_l, r0)*pmdec0*f_D*psi_pm
    J[..., 0, 5] = - pmra*(dt*f_T)**2/f_V - dot_last(p_l, r0)*f_D*psi_r

    # d(delta)/d(valores inicias)
    J[..., 1, 0] = dot_last(q_l, p0)*(1 + mi_r0*dt*f_T)*f_D - dot_last(q_l, r0)*pmra0*dt*f_T*f_D
    J[..., 1, 1] = dot_last(q_l, q0)*(1 + mi_r0*dt*f_T)*f_D - dot_last(q_l, r0)*pmdec0*dt*f_T*f_D
    J[..., 1, 2] = - dot_last(q_l, r0)*f_D*psi_parallax
    J[..., 1, 3] = dot_last(q_l, p0)*dt*f_T*f_D - dot_last(q_l, r0)*pmra0*f_D*psi_pm
    J[..., 1, 4] = dot_last(q_l, q0)*dt*f_T*f_D - dot_last(q_l, r0)*pmdec0*f_D*psi_pm
    J[..., 1, 5] = - pmdec*(dt*f_T)**2/f_V - dot_last(q_l, r0)*f_D*psi_r

    # d(parallax)/d(valores inicias)
    J[..., 2, 0] = 0
    J[..., 2, 1] = 0
    J[..., 2, 2] = f_D - parallax*(mi_r*dt*f_T/f_V)*psi_parallax
    J[..., 2, 3] = - parallax*pmra0*((dt*f_T)**2*f_D**2 + (mi_r*dt*f_T/f_V)*psi_pm)
    J[..., 2, 4] = - parallax*pmdec0*((dt*f_T)**2*f_D**2 + (mi_r*dt*f_T/f_V)*psi_pm)
    J[..., 2, 5] = - parallax*((1 + mi_r0*dt*f_T)*dt*f_T*f_D**2 + (mi_r*dt*f_T/f_V)*psi_r)

    # d(pmra)/d(valores inicias)
    J[..., 3, 0] = - dot_last(p_l, p0)*mi0**2*dt*f_T*f_D**3*f_V - dot_last(p_l, r0)*pmra0*(1+mi_r0*dt*f_T)*f_D**3*f_V
    J[..., 3, 1] = - dot_last(p_l, q0)*mi0**2*dt*f_T*f_D**3*f_V - dot_last(p_l, r0)*pmdec0*(1+mi_r0*dt*f_T)*f_D**3*f_V
    J[..., 3, 2] = dot_last(p_l, ni)*psi_parallax
    J[..., 3, 3] = dot_last(p_l, p0)*(1 + mi_r0*dt*f_T)*f_D**3*f_V - 2*dot_last(p_l, r0)*pmra0*dt*f_T*f_D**3*f_V \
        - 3*pmra*pmra0*(dt*f_T)**2*f_D**2*f_V + pmra*pmra0*chi_pm + dot_last(p_l, ni)*pmra0*psi_pm
    J[..., 3, 4] = dot_last(p_l, q0)*(1 + mi_r0*dt*f_T)*f_D**3*f_V - 2*dot_last(p_l, r0)*pmdec0*dt*f_T*f_D**3*f_V \
        - 3*pmra*pmdec0*(dt*f_T)**2*f_D**2*f_V + pmra*pmdec0*chi_pm + dot_last(p_l, ni)*pmdec0*psi_pm
    J[..., 3, 5] = dot_last(p_l, (vec_mi0*f_D[..., None] - 3*vec_mi*(1 + mi_r0*dt*f_T)[..., None]))*dt*f_T*f_D**2*f_V

    # d(pmdec)/d(valores inicias)
    J[..., 4, 0] = - dot_last(q_l, p0)*mi0**2*dt*f_T*f_D**3*f_V - dot_last(q_l, r0)*pmra0*(1+mi_r0*dt*f_T)*f_D**3*f_V
    J[..., 4, 1] = - dot_last(q_l, q0)*mi0**2*dt*f_T*f_D**3*f_V - dot_last(q_l, r0)*pmdec0*(1+mi_r0*dt*f_T)*f_D**3*f_V
    J[..., 4, 2] = dot_last(q_l, ni)*psi_parallax
    J[..., 4, 3] = dot_last(q_l, p0)*(1 + mi_r0*dt*f_T)*f_D**3*f_V - 2*dot_last(q_l, r0)*pmra0*dt*f_T*f_D**3*f_V \
        - 3*pmdec*pmra0*(dt*f_T)**2*f_D**2*f_V + pmdec*pmra0*chi_pm + dot_last(q_l, ni)*pmra0*psi_pm
    J[..., 4, 4] = dot_last(q_l, q0)*(1 + mi_r0*dt*f_T)*f_D**3*f_V - 2*dot_last(q_l, r0)*pmdec0*dt*f_T*f_D**3*f_V \
        - 3*pmdec*pmdec0*(dt*f_T)**2*f_D**2*f_V + pmdec*pmdec0*chi_pm + dot_last(q_l, ni)*pmdec0*psi_pm
    J[..., 4, 5] = dot_last(q_l, (vec_mi0*f_D[..., None] - 3*vec_mi*(1 + mi_r0*dt*f_T)[..., None]))*dt*f_T*f_D**2*f_V

    # d(rad_vel)/d(valores inicias)
    J[..., 5, 0] = 0
    J[..., 5, 1] = 0
    J[..., 5, 2] = eta*psi_parallax
    J[..., 5, 3] = 2*pmra0*(1 + mi_r0*dt*f_T)*dt*f_T*f_D**4*f_V + pmra0*mi_r*chi_pm + pmra0*eta*psi_pm
    J[..., 5, 4] = 2*pmdec0*(1 + mi_r0*dt*f_T)*dt*f_T*f_D**4*f_V + pmdec0*mi_r*chi_pm + pmdec0*eta*psi_pm
    J[..., 5, 5] = ((1 + mi_r0*dt*f_T)**2 - mi0**2*(dt*f_T)**2)*f_D**4*f_V + mi_r*chi_r + eta*psi_r

    # Propagated covariance matrix
    C = np.matmul(np.matmul(J, cov_matrix), np.swapaxes(J, -1, -2))

    err = np.sqrt(np.stack([C[..., i, i] for i in np.arange(3)], axis=-1))
    if not np.all(par_out):
        err = err[..., :2]

    return coord, err


def choice_star(catalogue, coord, columns, source):
    """Prompts the user to choose one star from a catalogue table.

    Parameters
    ----------
    catalogue : `astropy.table.Table`
        Catalogue table with candidate sources.
    coord : `astropy.coordinates.SkyCoord`
        Reference coordinate used to sort sources by distance.
    columns : `list`
        Column names to display. The first two columns must contain RA and DEC.
    source : `str`
        Source context used to decide how cancellation is handled.

    Returns
    -------
    catalogue : `astropy.table.Table`, `None`
        Table row selected by the user, or None when selection is cancelled for
        supported contexts.
    """
    from astropy.table import Table

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


def edr3ToICRF(pmra, pmdec, ra, dec, G):
    """Corrects Gaia EDR3 bright-star proper motions to the ICRF frame.

    Adapted from Cantat-Gaudin & Brandt, A&A, 2021.

    Parameters
    ----------
    pmra : `float`, `int`, `astropy.units.Quantity`
        Proper motion in right ascension, in mas/yr when unitless.

    pmdec : `float`, `int`, `astropy.units.Quantity`
        Proper motion in declination, in mas/yr when unitless.

    ra : `float`, `int`, `astropy.units.Quantity`
        Right ascension, in deg when unitless.

    dec : `float`, `int`, `astropy.units.Quantity`
        Declination, in deg when unitless.

    G : `float`, `int`
        Gaia G magnitude.

    Returns
    -------
    pmra_icrf, pmdec_icrf : `float`, `astropy.units.Quantity`
        Pair of corrected proper motions.
    """
    if G >= 13:
        return pmra, pmdec

    ra = u.Quantity(ra, unit=u.deg)
    dec = u.Quantity(dec, unit=u.deg)
    pmra = u.Quantity(pmra, unit=u.mas / u.year)
    pmdec = u.Quantity(pmdec, unit=u.mas / u.year)

    table1 = np.array([[0.0, 9.0, 18.4, 33.8, -11.3],
                       [9.0, 9.5, 14.0, 30.7, -19.4],
                       [9.5, 10.0, 12.8, 31.4, -11.8],
                       [10.0, 10.5, 13.6, 35.7, -10.5],
                       [10.5, 11.0, 16.2, 50.0, 2.1],
                       [11.0, 11.5, 19.4, 59.9, 0.2],
                       [11.5, 11.75, 21.8, 64.2, 1.0],
                       [11.75, 12.0, 17.7, 65.6, -1.9],
                       [12.0, 12.25, 21.3, 74.8, 2.1],
                       [12.25, 12.5, 25.7, 73.6, 1.0],
                       [12.5, 12.75, 27.3, 76.6, 0.5],
                       [12.75, 13.0, 34.9, 68.9, -2.9]]).T

    g_min = table1[0]
    g_max = table1[1]
    # pick the appropriate omegaXYZ for the source ’s magnitude :
    omega_x = table1[2][(g_min <= G) & (g_max > G)][0]*(u.mas/u.year)/1000.0
    omega_y = table1[3][(g_min <= G) & (g_max > G)][0]*(u.mas/u.year)/1000.0
    omega_z = table1[4][(g_min <= G) & (g_max > G)][0]*(u.mas/u.year)/1000.0
    pmra_corr = -1 * np.sin(dec) * np.cos(ra) * omega_x - np.sin(dec) * np.sin(ra) * omega_y + np.cos(dec) * omega_z
    pmdec_corr = np.sin(ra) * omega_x - np.cos(ra) * omega_y
    pmra_icrf  = pmra  - pmra_corr
    pmdec_icrf = pmdec - pmdec_corr
    return pmra_icrf, pmdec_icrf
