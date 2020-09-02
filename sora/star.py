from astropy.coordinates import SkyCoord, SphericalCosLatDifferential, Distance
from astropy.coordinates import get_sun, SphericalRepresentation, SkyOffsetFrame, ICRS
from astropy.coordinates import Longitude, Latitude
from astropy.time import Time
from astropy.table import Table
import astropy.units as u
import astropy.constants as const
from astroquery.vizier import Vizier
import warnings
import numpy as np
from sora.config import test_attr, input_tests


warnings.simplefilter('always', UserWarning)


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
        print('Downloading star parameters from {}'.format(kwargs['catalog']))
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
    return catalogue[[k[choice-1]]]


class Star():
    def __init__(self, **kwargs):
        """ Defines a star

        Parameters:
            code (str): Gaia-DR2 Source code for searching in VizieR.
            coord (str, SkyCoord): if code is not given, coord nust have the coordinates
                RA and DEC of the star to search in VizieR: 'hh mm ss.ss +dd mm ss.ss'
            ra (int, float): Right Ascension, in deg.
            dec (int, float): Declination, in deg.
            parallax (int, float): Parallax, in mas. Default = 0
            pmra (int, float): Proper Motion in RA*, in mas/year. Default = 0
            pmdec (int, float): Proper Motion in DEC, in mas/year. Default = 0
            rad_vel (int, float): Radial Velocity, in km/s. Default = 0 km/s
            epoch (str, Time): Epoch of the coordinates. Default = 'J2000'
            nomad (bool): If true, it tries to download the magnitudes from NOMAD
                catalogue. Default = True
            log (bool): If true, it prints the downloaded information. Default = True
            local (bool): If true, it uses the given coordinate in 'coord'
                as final coordinate. Default = False

        - The user can only give ("coord") or ("ra" and "dec"), but not both.
        - To download the coordinates from Gaia, "local" must be set as False
            and the ("code") or ("coord") or ("ra" and "dec") must be given.
            All values downloaded from Gaia will replace the ones given by the user.
        """
        self.__attributes = {}
        self.mag = {}
        self.errors = {'RA': 0, 'DEC': 0, 'Plx': 0, 'pmRA': 0, 'pmDE': 0, 'rad_vel': 0}
        allowed_kwargs = ['code', 'coord', 'dec', 'epoch', 'local', 'log', 'nomad', 'parallax', 'pmdec', 'pmra', 'ra', 'rad_vel']
        input_tests.check_kwargs(kwargs, allowed_kwargs=allowed_kwargs)
        self.__log = kwargs.get('log', True)
        local = kwargs.get('local', False)
        if 'code' in kwargs:
            self.code = kwargs['code']
        if 'coord' in kwargs:
            if 'ra' in kwargs or 'dec' in kwargs:
                raise ValueError("User must give 'coord' or 'ra' and 'dec', not both")
            coord = SkyCoord(kwargs['coord'], unit=('hourangle', 'deg'))
            self.ra = coord.ra
            self.dec = coord.dec
        if 'ra' in kwargs and 'dec' in kwargs:
            self.ra = kwargs.get('ra')
            self.dec = kwargs.get('dec')
        self.parallax = kwargs.get('parallax', 0.0)
        self.pmra = kwargs.get('pmra', 0.0)
        self.pmdec = kwargs.get('pmdec', 0.0)
        self.rad_vel = kwargs.get('rad_vel', 0.0)
        self.epoch = kwargs.get('epoch', 'J2000')
        if local:
            if 'RA' not in self.__attributes or 'DEC' not in self.__attributes:
                raise ValueError("User must give 'ra' and 'dec' for local coordinates")
        else:
            if not hasattr(self, 'code') and 'RA' not in self.__attributes:
                raise ValueError("User must give gaia Source ID 'code' or coordinates for the online search")
            self.__searchgaia()
        if kwargs.get('nomad', True):
            self.__getcolors()

    @property
    def ra(self):
        """
        Return the Right Ascension of the star
        """
        return self.__attributes['RA']

    @ra.setter
    def ra(self, value):
        """ Define Right Ascension

        Parameter:
            value(int, float, astropy.unit): RA in deg
        """
        self.__attributes['RA'] = Longitude(value, unit=u.hourangle)

    @property
    def dec(self):
        """
        Return the Declination of the star
        """
        return self.__attributes['DEC']

    @dec.setter
    def dec(self, value):
        """ Define Declination

        Parameter:
            value(int, float, astropy.unit): DEC in deg
        """
        self.__attributes['DEC'] = Latitude(value, unit=u.deg)

    @property
    def parallax(self):
        """
        Return the Parallax of the star
        """
        return self.__attributes.get('PAR', 0*u.mas)

    @parallax.setter
    def parallax(self, value):
        """ Define Parallax

        Parameter:
            value(int, float, astropy.unit): Parallax in mas
        """
        par = u.Quantity(value, unit=u.mas)
        if par <= 0*u.mas:
            par = 0*u.mas
        self.__attributes['PAR'] = par

    @property
    def distance(self):
        """
        Return the Distance of the star
        """
        if self.parallax > 0*u.mas:
            return Distance(parallax=self.parallax, allow_negative=False)
        else:
            raise ValueError('SORA is not able to determine distance from paralax {}'.format(self.parallax))

    @property
    def coord(self):
        try:
            return SkyCoord(self.ra, self.dec, self.distance)
        except ValueError:
            return SkyCoord(self.ra, self.dec)

    @property
    def pmra(self):
        """
        Return the Proper Motion in Right Ascension of the star
        """
        return self.__attributes.get('PMRA', 0*u.mas/u.year)

    @pmra.setter
    def pmra(self, value):
        """ Define Parallax

        Parameter:
            value(int, float, astropy.unit): RA Proper Motion in mas/year
        """
        self.__attributes['PMRA'] = u.Quantity(value, unit=u.mas/u.year)

    @property
    def pmdec(self):
        """
        Return the Proper Motion in Declination of the star
        """
        return self.__attributes.get('PMDEC', 0*u.mas/u.year)

    @pmdec.setter
    def pmdec(self, value):
        """ Define Parallax

        Parameter:
            value(int, float, astropy.unit): DEC Proper Motion in mas/year
        """
        self.__attributes['PMDEC'] = u.Quantity(value, unit=u.mas/u.year)

    @property
    def rad_vel(self):
        """
        Return the Radial Velocity of the star
        """
        return self.__attributes.get('RAD_VEL', 0*u.km/u.s)

    @rad_vel.setter
    def rad_vel(self, value):
        """ Define Radial Velocity

        Parameter:
            value(int, float, astropy.unit): Radial Velocity in km/s
        """
        self.__attributes['RAD_VEL'] = u.Quantity(value, unit=u.km/u.s)

    @property
    def epoch(self):
        """
        Return the Epoch of the position of the star
        """
        return self.__attributes['EPOCH']

    @epoch.setter
    def epoch(self, value):
        """ Define Radial Velocity

        Parameter:
            value(int, float, astropy.unit): Radial Velocity in km/s
        """
        self.__attributes['EPOCH'] = Time(value)

    def set_magnitude(self, **kwargs):
        """ Sets the magnitudes of a star.

        Parameters:
            (band name)=(float): The star magnitude for given band. The band name can be any string the user wants.

        Examples:
            set_magnitude(G=10)
            set_magnitude(K=15)
            set_magnitude(newband=6)
        """
        for key in kwargs:
            mag = test_attr(kwargs[key], float, key)
            if key in self.mag:
                warnings.warn('{0} mag already defined. {0}={1} will be replaced by {0}={2}'.format(
                    key, self.mag[key], mag))
            self.mag[key] = mag

    def set_diameter(self, diameter):
        """ Sets an user diameter for the star, in mas.

        Parameters:
            diameter (int,float): sets the user diameter of the star, in mas
        """
        self.diameter_user = diameter*u.mas

    def van_belle(self):
        """ Determines the diameter of a star in mas using equations from van Belle (1999)
            -- Publi. Astron. Soc. Pacific 111, 1515-1523:
        """
        return van_belle(self.mag.get('B'), self.mag.get('V'), self.mag.get('K'))

    def kervella(self):
        """ Determines the diameter of a star in mas using equations from Kervella et. al (2004)
            -- A&A Vol.  426, No.  1:
        """
        return kervella(self.mag.get('B'), self.mag.get('V'), self.mag.get('K'))

    def apparent_diameter(self, distance, mode='auto', band='V', star_type='sg', log=True):
        """ Calculates the apparent diameter of the star at a given distance

        Parameters:
            distance (int, float): Object geocentric distance, in AU
            mode (str): The mode to calculate the apparent diameter
                'user': calculates using user given diameter
                'gaia': calculates using diameter obtained from Gaia
                'kervella': calculates using Kervella equations
                'van_belle': calculates using van Belle equations
                'auto' (default): tries all the above methods until it is able to calculate diameter.
                    The order of try is the same as shown above (user, Gaia, Kervella, Van Belle).
            'band' (str): The band filter to calculate the diameter.
                If mode is 'kervella' or 'van_belle', the filter must be given, 'B' or 'V'.
                If mode 'auto', 'V' is selected.
            'star_type' (str): type of star to calculate the diameter.
                If mode is 'van_belle', the star type must be given.
                If mode is 'auto', type = 'sg'.
                Types can be:
                    - 'sg' for 'Super Giant'
                    - 'ms' for 'Main Sequence'
                    - 'vs' for 'Variable Star'
            'log' (bool): If True, it prints the mode used by 'auto'.
        """
        try:
            distance = distance.to(u.km)
        except:
            distance = distance*u.AU

        if mode in ['user', 'auto']:
            try:
                diam = distance*np.tan(self.diameter_user)
                if log:
                    print('Calculating apparent diameter from user defined diameter')
                return diam.to(u.km)
            except:
                pass

        if mode == 'user':
            raise ValueError('User diameter must be informed.')

        if mode in ['gaia', 'auto']:
            try:
                diam = distance*np.tan(self.diameter_gaia)
                if log:
                    print('Apparent diameter using Gaia')
                return diam.to(u.km)
            except:
                pass

        if mode == 'gaia':
            raise ValueError('It is not possible to calculate star diameter from Gaia.')

        if band not in ['B', 'V']:
            raise KeyError('band must be informed as "B", or "V"')

        if mode in ['kervella', 'auto']:
            diam_kerv = self.kervella().get(band)
            if diam_kerv is None:
                raise ValueError('Diameter could not be calculated for given band')
            if log:
                print('Apparent diameter using Kervella et al. (2004)')
            diam = distance*np.tan(diam_kerv)
            return diam.to(u.km)

        if star_type not in ['sg', 'ms', 'vs']:
            raise KeyError('star_type must be informed as "sg", "ms" or "vs"')

        if mode in ['van_belle', 'auto']:
            diam_van = self.van_belle().get(star_type)
            if diam_van is None:
                raise ValueError('Diameter could not be calculated using Van Belle')
            diam_van = diam_van.get(band)
            if diam_van is None:
                raise ValueError('Diameter could not be calculated for given band')
            if log:
                print('Apparent diameter using van Belle (1999)')
            diam = distance*np.tan(diam_van)
            return diam.to(u.km)

        raise AttributeError("Star apparent diameter could not be calculated. ",
                             "Please define star diameter or B,V,K magnitudes.")

    def __searchgaia(self):
        """ Searches for the star position in the Gaia catalogue and save informations
        """
        if hasattr(self, 'code'):
            catalogue = search_star(code=self.code, columns=['**'], catalog='I/345/gaia2', log=self.__log)
        else:
            catalogue = search_star(coord=self.coord, columns=['**'], radius=2*u.arcsec,
                                    catalog='I/345/gaia2', log=self.__log)
        if len(catalogue) == 0:
            raise ValueError('No star was found. It does not exist or VizieR is out.')
        catalogue = catalogue[0]
        if len(catalogue) > 1:
            if self.__log:
                print('{} stars were found within 2 arcsec from given coordinate.'.format(len(catalogue)))
                print('The list below is sorted by distance. Please select the correct star')
            catalogue = choice_star(catalogue, self.coord, ['RA_ICRS', 'DE_ICRS', 'Gmag'])
        self.code = catalogue['Source'][0]
        self.ra = catalogue['RA_ICRS'][0]*u.deg
        self.dec = catalogue['DE_ICRS'][0]*u.deg
        self.pmra = catalogue['pmRA'][0]*u.mas/u.year
        self.pmdec = catalogue['pmDE'][0]*u.mas/u.year
        self.epoch = Time(catalogue['Epoch'][0], format='jyear')
        self.parallax = catalogue['Plx'][0]*u.mas
        self.rad_vel = catalogue['RV'][0]*u.km/u.s
        self.set_magnitude(G=catalogue['Gmag'][0])

        self.meta_gaia = {c: catalogue[c][0] for c in catalogue.columns}
        rad = catalogue['Rad'][0]
        if np.ma.core.is_masked(rad) or np.ma.core.is_masked(catalogue['Plx'][0]):
            if self.__log:
                warnings.warn('Gaia catalogue does not have the star radius.')
        else:
            self.diameter_gaia = 2*np.arctan((rad*u.solRad)/self.distance).to(u.mas)

        self.errors['RA'] = self.meta_gaia['e_RA_ICRS']*u.mas
        self.errors['DEC'] = self.meta_gaia['e_DE_ICRS']*u.mas
        self.errors['Plx'] = self.meta_gaia['e_Plx']*u.mas
        self.errors['pmRA'] = self.meta_gaia['e_pmRA']*(u.mas/u.yr)
        self.errors['pmDE'] = self.meta_gaia['e_pmDE']*(u.mas/u.yr)
        self.errors['rad_vel'] = self.meta_gaia['e_RV']*(u.km/u.s)

        A = (1*u.AU).to(u.km).value
        cov = np.zeros((6, 6))
        a = ['RA', 'DE', 'Plx', 'pmRA', 'pmDE']
        for i in np.arange(5):
            v1 = 'e_' + a[i]
            if i in [0, 1]:
                v1 += '_ICRS'
            for j in np.arange(i, 5):
                v2 = 'e_' + a[j]
                if j in [0, 1]:
                    v2 += '_ICRS'
                if i == j:
                    x = self.meta_gaia[v1]**2
                    if not np.ma.core.is_masked(x):
                        cov[i, i] = x
                else:
                    x = self.meta_gaia[a[i]+a[j]+'cor']*self.meta_gaia[v1]*self.meta_gaia[v2]
                    if not np.ma.core.is_masked(x):
                        cov[i, j] = x
                        cov[j, i] = cov[i, j]
            x = cov[i, 2]*(self.meta_gaia['RV']/A)
            if not np.ma.core.is_masked(x):
                cov[i, 5] = x
                cov[5, i] = cov[i, 5]
        x = cov[2, 2]*(self.meta_gaia['RV']**2 + self.meta_gaia['e_RV']**2)/(A**2) \
            + (self.meta_gaia['Plx']*self.meta_gaia['e_RV']/A)**2
        if not np.ma.core.is_masked(x):
            cov[5, 5] = x
        cov[np.where(np.isnan(cov))] = 0.0
        self.cov = cov

        if self.__log:
            print('1 Gaia-DR2 star found G={}'.format(catalogue['Gmag'][0]))
            print('star coordinate at J{}: RA={} +/- {}, DEC={} +/- {}'.format(self.epoch.jyear,
                  self.ra.to_string(u.hourangle, sep='hms', precision=5), self.errors['RA'],
                  self.dec.to_string(u.deg, sep='dms', precision=4), self.errors['DEC']))

    def __getcolors(self):
        """ Searches for the B,V,K magnitudes of the star in the NOMAD catalogue on VizieR
        """
        columns = ['RAJ2000', 'DEJ2000', 'Bmag', 'Vmag', 'Rmag', 'Jmag', 'Hmag', 'Kmag']
        catalogue = search_star(coord=self.coord, columns=columns, radius=2*u.arcsec,
                                catalog='I/297/out', log=self.__log)
        if len(catalogue) == 0:
            if self.__log:
                warnings.warn('No star was found on NOMAD that matches the star')
            return
        catalogue = catalogue[0]
        if len(catalogue) > 1:
            print('{} stars were found within 2 arcsec from given coordinate.'.format(len(catalogue)))
            print('The list below is sorted by distance. Please select the correct star')
            if hasattr(self.mag, 'G'):
                print('Star G mag: {}'.format(self.mag['G']))
            catalogue = choice_star(catalogue, self.coord, ['RAJ2000', 'DEJ2000', 'Bmag', 'Vmag',
                                                            'Rmag', 'Jmag', 'Hmag', 'Kmag'])
            if catalogue is None:
                return
        errors = []
        for mag in ['B', 'V', 'R', 'J', 'H', 'K']:
            name = mag + 'mag'
            if np.ma.core.is_masked(catalogue[name][0]):
                errors.append(mag)
                continue
            self.set_magnitude(**{mag: catalogue[name][0]})
        if len(errors) > 0 and self.__log:
            print('Magnitudes in {} were not located in NOMAD'.format(errors))

    def geocentric(self, time):
        """ Calculates the position of the star, propagating the position using parallax and proper motion

        Parameters:
            time (float, Time): reference time to apply proper motion and calculate paralax.
        """
        try:
            time = Time(time)
        except:
            time = Time(time, format='jd', scale='utc')
        n_coord = self.barycentric(time)
        if self.coord.distance.unit.is_unity() or np.isnan(self.coord.distance):
            g_coord = n_coord
        else:
            sun = get_sun(time)
            g_coord = SkyCoord(*(n_coord.cartesian.xyz + sun.cartesian.xyz),
                               representation_type='cartesian')
            g_coord = g_coord.represent_as(SphericalRepresentation)
            g_coord = SkyCoord(g_coord.lon, g_coord.lat, g_coord.distance)

        if hasattr(self, 'offset'):
            star_frame = SkyOffsetFrame(origin=g_coord)
            new_pos = SkyCoord(lon=self.offset.d_lon_coslat, lat=self.offset.d_lat, frame=star_frame)
            return new_pos.transform_to(ICRS)

        return g_coord

    def barycentric(self, time):
        """ Calculates the position of the star using proper motion

        Parameters:
            time (str, Time): reference time to apply proper motion.
        """
        try:
            time = Time(time)
        except:
            time = Time(time, format='jd', scale='utc')
        dt = time - self.epoch
        n_coord = spatial_motion(self.ra, self.dec, self.pmra, self.pmdec, self.parallax, self.rad_vel,  dt=dt.jd)
        return n_coord

    def error_at(self, time):
        """ Estimates the star position error at a given time

        Parameters:
            time (str, Time): reference time to project star error.

        Returns:
            errors in RA* and DEC
        """
        try:
            time = Time(time)
        except:
            time = Time(time, format='jd', scale='utc')
        dt = time - self.epoch
        n_coord, errors = spatial_motion(self.ra, self.dec, self.pmra, self.pmdec, self.parallax,
                                         self.rad_vel,  dt=dt.jd, cov_matrix=self.cov)
        return errors[0]*u.mas, errors[1]*u.mas

    def add_offset(self, da_cosdec, ddec):
        """ Adds an offset to the star position

        Parameters:
            da_cosdec (int, float): offset in Delta_alpha_cos_delta in mas
            ddec (int, float): offset in Delta_delta in mas
        """
        dadc = test_attr(da_cosdec, float, 'da_cosdec')
        dd = test_attr(ddec, float, 'ddec')
        self.offset = SphericalCosLatDifferential(dadc*u.mas, dd*u.mas, 0.0*u.km)

    def __str__(self):
        """ String representation of the Star class
        """
        out = ''
        if hasattr(self, 'code'):
            out += 'Gaia-DR2 star Source ID: {}\n'.format(self.code)
        else:
            out += 'User coordinates\n'
        out += ('ICRS star coordinate at J{}:\n'
                'RA={} +/- {:.4f}, DEC={} +/- {:.4f}\n'
                'pmRA={:.3f} +/- {:.3f} mas/yr, pmDEC={:.3f} +/- {:.3f} mas/yr\n'
                'Plx={:.4f} +/- {:.4f} mas, Rad. Vel.={:.2f} +/- {:.2f} km/s \n\n'.format(
                    self.epoch.jyear, self.ra.to_string(u.hourangle, sep='hms', precision=5),
                    self.errors['RA'], self.dec.to_string(u.deg, sep='dms', precision=4), self.errors['DEC'],
                    self.pmra.value, self.errors['pmRA'], self.pmdec.value, self.errors['pmDE'],
                    self.parallax.value, self.errors['Plx'], self.rad_vel.value, self.errors['rad_vel']))
        if hasattr(self, 'offset'):
            out += 'Offset Apllied: d_alpha_cos_dec = {}, d_dec = {}\n'.format(
                self.offset.d_lon_coslat, self.offset.d_lat)
        out += 'Magnitudes:'
        mag_out = [' {}: {:6.3f}'.format(mag, self.mag[mag]) for mag in self.mag]
        out_mag = []
        for i, mag in enumerate(mag_out):
            if i % 6 == 0:
                out_mag.append([])
            out_mag[-1].append(mag)
        out += (',\n'+' '*11).join([','.join(out_i) for out_i in out_mag])
        out += '\n\n'
        if hasattr(self, 'diameter_gaia'):
            out += 'Apparent diameter: {:.4f}, Source: Gaia-DR2\n'.format(self.diameter_gaia)
        if hasattr(self, 'diameter_user'):
            out += 'Apparent diameter: {:.4f}, Source: User\n'.format(self.diameter_user)
        kerv = self.kervella()
        if kerv:
            out += 'Apparent diameter from Kervella et. al (2004):\n'
            out += '   ' + ','.join([' {}: {:.4f}'.format(k, v) for k, v in kerv.items()])
        vanb = self.van_belle()
        if vanb:
            out += '\nApparent diameter from van Belle (1999):'
            for key, value in vanb.items():
                out += '\n    {}:'.format(key)
                out += ','.join([' {}: {:.4f}'.format(k, v) for k, v in value.items()])
        return out
