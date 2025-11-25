import astropy.units as u
import numpy as np

__all__ = ['project_to_ring_plane','to_sky_plane','calc_longitude','calc_equatorial_vel', 'calc_coef_projecao']

def calc_coef_projecao(ephem, pole_coord, B, P, earth_pole):
    """ Function that calculate the projection coeficients """ 

    east    = np.array([-np.sin(ephem.ra), np.cos(ephem.ra), 0])
    north   = np.array([-np.cos(ephem.ra)*np.sin(ephem.dec), -np.sin(ephem.ra)*np.sin(ephem.dec), np.cos(ephem.dec)])
    target  = np.array([np.cos(ephem.ra)*np.cos(ephem.dec), np.sin(ephem.ra)*np.cos(ephem.dec), np.sin(ephem.dec)])
    pt = np.array([np.cos(earth_pole.ra)*np.cos(earth_pole.dec), np.sin(earth_pole.ra)*np.cos(earth_pole.dec), np.sin(earth_pole.dec)])

    east_pole = np.cos(pole_coord.dec)*np.sin(pole_coord.ra - ephem.ra)
    north_pole = np.sin(pole_coord.dec)*np.cos(ephem.dec) - np.cos(pole_coord.dec)*np.sin(ephem.dec)*np.cos(pole_coord.ra - ephem.ra)

    a = -east*np.sin(B)*np.sin(P) - north*np.sin(B)*np.cos(P) - target*np.cos(B)
    b = -east*np.cos(P) + north*np.sin(P)
    p = +east*np.cos(B)*np.sin(P) + north*np.cos(B)*np.cos(P) - target*np.sin(B)

    n1 = np.cross(pt, p)
    module_n1 = np.sqrt(np.sum(n1**2))
    n1 = n1/module_n1
    n2 = np.cross(p, n1)

    coef_polo = np.array([np.dot(a,n1), np.dot(b,n1), np.dot(a,n2), np.dot(b,n2)])
    coef = np.array([
        -east_pole/(np.sin(B)*np.cos(B)),
        -north_pole/(np.sin(B)*np.cos(B)),
        -north_pole/np.cos(B),
        +east_pole/np.cos(B)
    ])

    return coef, coef_polo


def project_to_ring_plane(ksi, eta, coef, coef_polo, ksi_0=0.0, eta_0=0.0):
    """ Convert between the sky-plane and the (equatorial) ring plane 

        Parameters:
            ksi (int, float): Ksi
            eta (int, float): Eta
            ksi_0 (int, float): Central ksi value, default=0.0
            eta_0 (int, float): Central eta value, default=0.0
    
        Returns:
            xp (float): X Position in the ring plane
            yp (float): Y Position in the ring plane            
        """
    delta_ksi = ksi - ksi_0
    delta_eta = eta - eta_0

    x = delta_ksi*coef[0] + delta_eta*coef[1]
    y = delta_ksi*coef[2] + delta_eta*coef[3]
    xp = x*coef_polo[0] + y*coef_polo[1]
    yp = x*coef_polo[2] + y*coef_polo[3]
    return xp, yp


def to_sky_plane(xp, yp, coef, coef_polo, ksi_0=0.0, eta_0=0.0):
    """ Convert between the sky-plane and the (equatorial) ring plane 

        Parameters:
            xp (int, float): X Position in the ring plane
            yp (int, float): Y Position in the ring plane
            ksi_0 (int, float): Central ksi value, default=0.0
            eta_0 (int, float): Central eta value, default=0.0

        Returns:
            ksi (float): Ksi
            eta (float): Eta            
        """
    denom = coef_polo[0]*coef_polo[3] - coef_polo[2]*coef_polo[1]
    x = (coef_polo[3]*xp - coef_polo[1]*yp) / denom
    y = (yp / coef_polo[3]) - (x * coef_polo[2] / coef_polo[3])

    denom2 = coef[0]*coef[3] - coef[2]*coef[1]
    delta_ksi = (coef[3]*x - coef[1]*y) / denom2
    delta_eta = (y / coef[3]) - (delta_ksi * coef[2] / coef[3])

    ksi = delta_ksi + ksi_0
    eta = delta_eta + eta_0
    return ksi, eta


def calc_longitude(ksi, eta, coef, coef_polo, ksi_0=0.0, eta_0=0.0):
    """ Calculate the True longitude in the J2000 reference frame 

        Parameters:
            ksi (int, float): Ksi
            eta (int, float): Eta
            ksi_0 (int, float): Central ksi value, default=0.0
            eta_0 (int, float): Central eta value, default=0.0

        Returns:
            along (float): True longitude in the J2000 reference frame, in degrees
        """
    delta_ksi = ksi - ksi_0
    delta_eta = eta - eta_0

    x = delta_ksi*coef[0] + delta_eta*coef[1]
    y = delta_ksi*coef[2] + delta_eta*coef[3]
    xp = x*coef_polo[0] + y*coef_polo[1]
    yp = x*coef_polo[2] + y*coef_polo[3]

    along = np.arctan2(yp, xp) * u.rad
    return along.to('deg').value


def calc_equatorial_vel(ksi, eta, vel_ksi, vel_eta, coef, coef_polo, ksi_0=0.0, eta_0=0.0):
    """ Calculate the radial velocity in the equatorial (ring) plane

        Parameters:
            ksi (int, float): Ksi
            eta (int, float): Eta
            vel_ksi (int, float): Velocity Ksi
            vel_eta (int, float): Velocity Eta
            ksi_0 (int, float): Central ksi value, default=0.0
            eta_0 (int, float): Central eta value, default=0.0

        Returns:
            v_rp (float): Radial velocity in the equatorial (ring) plane, in km/s
        """
    delta_ksi = ksi - ksi_0
    delta_eta = eta - eta_0

    x = delta_ksi*coef[0] + delta_eta*coef[1]
    y = delta_ksi*coef[2] + delta_eta*coef[3]
    xp = x*coef_polo[0] + y*coef_polo[1]
    yp = x*coef_polo[2] + y*coef_polo[3]
    rp = np.sqrt(xp**2 + yp**2)

    xf = (delta_ksi + vel_ksi)*coef[0] + (delta_eta + vel_eta)*coef[1]
    yf = (delta_ksi + vel_ksi)*coef[2] + (delta_eta + vel_eta)*coef[3]
    xpf = xf*coef_polo[0] + yf*coef_polo[1]
    ypf = xf*coef_polo[2] + yf*coef_polo[3]
    rpf = np.sqrt(xpf**2 + ypf**2)

    v_rp = rpf - rp
    return v_rp
