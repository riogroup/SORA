from astropy.coordinates import SkyCoord, SphericalCosLatDifferential, Distance
from astropy.time import Time
import astropy.units as u
from astroquery.vizier import Vizier
import warnings
import numpy as np
from .config import test_attr

### Object for star
def search_star(**kwargs):
    """ Search position on Vizier and returns a catalogue.

    Parameters:
        coord (str, astropy.SkyCoord):Coordinate to search around.
        code (str):Gaia Source_id of the star
        columns (list):List of strings with the name of the columns to retrieve.
        radius (int, float, unit.quantity):Radius to search around coord.
        catalog (str):Vizier catalogue to search.

    Returns:
        catalogue(astropy.Table):An astropy Table with the catalogue informations.   

    """
    row_limit = 100
    print('Downloading star parameters from {}'.format(kwargs['catalog']))
    vquery = Vizier(columns=kwargs['columns'], row_limit=row_limit, timeout=600)
    if 'code' in kwargs:
        catalogue = vquery.query_constraints(catalog=kwargs['catalog'], Source=kwargs['code'])
    elif 'coord' in kwargs:
        catalogue = vquery.query_region(kwargs['coord'], radius=kwargs['radius'], catalog=kwargs['catalog'])
    return catalogue

def van_belle(magB, magV, magK):
    '''
    Determine the diameter of a star in mas using equations from van Belle (1999) 
    -- Publi. Astron. Soc. Pacific 111, 1515-1523:
    Inputs:
    magB, magV, magK: The magnitudes B, V and K of the star
    '''

    def calc_diameter(a1, a2, mag):
        return 10**(a1 + a2*(mag - magK) - 0.2*mag)

    params = {'sg': {'B': [ 0.648, 0.220], 'V': [0.669, 0.223]},
              'ms': {'B': [ 0.500, 0.290], 'V': [0.500, 0.264]},
              'vs': {'B': [ 0.789, 0.218], 'V': [0.840, 0.211]}}

    mag = np.array([magB, magV])
    diameter = {}
    for st in ['sg', 'ms', 'vs']:
        diameter[st] = {}
        for i,m in enumerate(['B','V']):
            diameter[st][m] = calc_diameter(*params[st][m], mag[i])*u.mas
    return diameter
    
def kervella(magB, magV, magK):
    '''
    Determine the diameter of a star in mas using equations from Kervella et. al (2004) 
    -- A&A Vol.  426, No.  1:
    Inputs:
    magB, magV, magK: The magnitudes B, V and K of the star
    '''
    const1 = np.array([0.0755, 0.0535])
    const2 = np.array([0.5170, 0.5159])
    mag = np.array([magV,magB])
    vals = 10**(const1*(mag-magK)+const2-0.2*magK)
    return {'V': vals[0]*u.mas, 'B': vals[1]*u.mas}

class Star():
    def __init__(self,**kwargs):
        '''
        Docstring
        Define a star
        '''
        self.__local = False
        self.mag = {}
        self.errors = {}
        if 'local' in kwargs:
            self.__local = test_attr(kwargs['local'], bool, 'local')
        if not any(i in kwargs for i in ['coord', 'code']):
            raise KeyError('Input values must have either the coordinate of the star or the Gaia code for online search')
        if 'coord' in kwargs:
            self.coord = SkyCoord(kwargs['coord'], unit=('hourangle', 'deg'))
            #self.coord = test_attr(kwargs['coord'], SkyCoord, 'coord')
        else:
            self.__local=False
        if 'code' in kwargs:
            self.code = test_attr(kwargs['code'], str, 'code')
        if not self.__local:
            self.__searchgaia2()
    
    
    def set_magnitude(self,**kwargs):
        '''
        Set the magnitudes of a star in the G, B, V and K band.
        usually this values can be found in the GDR2 and in the NOMAD catalogue.
        Inputs:
        radius = float, in mas
        '''
        for key in kwargs:
            mag = test_attr(kwargs[key], float, key)
            if key in self.mag:
                warnings.warn('{0} mag already defined. {1} will be replaced by {2}'.format(key, self.mag[key], kwargs[key]))
            self.mag[key] = kwargs[key]

            
    def set_diameter(self, diameter):
        '''
        Set the diameter of a star in mas.
        Inputs:
        star_diameter = float, in mas
        '''
        self.diameter_user = star_diameter*u.mas

    
    def van_belle(self):
        '''
        Determine the diameter of a star in mas using equations from van Belle (1999) 
        -- Publi. Astron. Soc. Pacific 111, 1515-1523:
        Inputs:
        magB, magV, magK: The magnitudes B, V and K of the star
        '''
        return van_belle(self.mag['B'], self.mag['V'], self.mag['K'])
    
    
    def kervella(self):
        '''
        Determine the diameter of a star in mas using equations from Kervella et. al (2004) 
        -- A&A Vol.  426, No.  1:
        Inputs:
        magB, magV, magK: The magnitudes B, V and K of the star
        '''
        return kervella(self.mag['B'], self.mag['V'], self.mag['K'])
    
    
    def apparent_diameter(self, distance, mode='auto', **kwargs):
        # calculate the apparent radius of the star at given distance
        if mode == 'auto':
            kwargs['obs_filter'] = 'V'
            kwargs['star_type'] = 'sg'
            
        if mode in ['user', 'auto']:
            try:
                diam = distance*np.tan(self.diameter_user)
                print('Calculating apparent diameter from user defined diameter')
                return diam.to(u.km)
            except:
                pass
        
        if mode == 'user':
            raise ValueError('User diameter must be informed.')
            
        if mode in ['gaia', 'auto']:
            try:
                diam = distance*np.tan(self.diameter_gaia)
                print('Apparent diameter using Gaia')
                return diam.to(u.km)
            except:
                pass
                
        if mode == 'gaia':
            raise ValueError('It is not possible to calculate star diameter from Gaia.')
                
        if 'obs_filter' not in kwargs:
            raise KeyError('obs_filter must be informed as "B", or "V"')
            
        if mode in ['kervella', 'auto']:
            if all(i in self.mag for i in ['B', 'V', 'K']):
                print('Apparent diameter using Kervella et al. (2004)')
                diam = distance*np.tan(self.kervella()[kwargs['obs_filter']])
                return diam.to(u.km)
            
        if 'star_type' not in kwargs:
            raise KeyError('star_type must be informed as "sg", "ms" or "vs"')
            
        if mode in ['van_belle', 'auto']:
            if all(i in self.mag for i in ['B', 'V', 'K']):
                print('Apparent diameter using van Belle (1999)')
                diam = distance*np.tan(self.van_belle()[kwargs['star_type']][kwargs['obs_filter']])
                return diam.to(u.km)
                        
        raise AttributeError('Star apparent diameter could not be calculated.\
Please define star diameter or B,V,K magnitudes.')
        
    
    def __searchgaia2(self):
        """search for the star position in the gaia catalogue and save informations
        """
        columns = ['Source', 'RA_ICRS', 'e_RA_ICRS', 'DE_ICRS', 'e_DE_ICRS', 'Plx', 'pmRA', 'e_pmRA', 'pmDE', 'e_pmDE', 'Gmag', 'e_Gmag', 'Dup', 'Epoch', 'Rad']
        if hasattr(self, 'code'):
            catalogue = search_star(code=self.code, columns=columns, catalog='I/345/gaia2')[0]
        else:
            catalogue = search_star(coord=self.coord, columns=columns, radius=2*u.arcsec, catalog='I/345/gaia2')[0]
        if len(catalogue) == 1:
            self.code = catalogue['Source']
            ra = catalogue['RA_ICRS']
            dec = catalogue['DE_ICRS']
            distance = Distance(parallax=catalogue['Plx'].quantity, allow_negative=True)
            pmra = catalogue['pmRA']
            pmde = catalogue['pmDE']
            epoch = Time(catalogue['Epoch'].quantity, format='jyear')
            self.coord = SkyCoord(ra, dec, distance=distance, pm_ra_cosdec=pmra, pm_dec=pmde, obstime=epoch)[0]
            self.set_magnitude(G=catalogue['Gmag'][0])
            self.errors['RA'] = catalogue['e_RA_ICRS'][0]*u.mas
            self.errors['DEC'] = catalogue['e_DE_ICRS'][0]*u.mas
            self.errors['pmRA'] = catalogue['e_pmRA'][0]*(u.mas/u.yr)
            self.errors['pmDEC'] = catalogue['e_pmDE'][0]*(u.mas/u.yr)
            rad = catalogue['Rad'][0]
            if np.ma.core.is_masked(rad):
                warnings.warn('Gaia catalogue does not have star radius.')
            else:
                self.diameter_gaia = 2*np.arctan((rad*u.solRad)/distance[0]).to(u.mas)
            self.__getcolors()
        else:
            ## pegar todas as estrelas e colocar como opcao pro usuario escolher.
            warnings.warn('{} stars found in the region searched'.format(len(catalogue)))
        
        
    def __getcolors(self):
        # search for the B,V,K magnitudes of the star on Vizier and saves the result
        columns = ['RAJ2000', 'DEJ2000', 'Bmag', 'Vmag', 'Rmag', 'Jmag', 'Hmag', 'Kmag']
        catalogue = search_star(coord=self.coord, columns=columns, radius=2*u.arcsec, catalog='I/297/out')[0]
        if len(catalogue) == 0:
            raise IndexError('No star was found on NOMAD that matches the star')
        elif len(catalogue) > 1:
            print('One or more stars were found in the region. Please select the correct one')
            n_coord = SkyCoord(catalogue['RAJ2000'], catalogue['DEJ2000'])
            ### pegar todas as estrelas e colocar como opcao pro usuario escolher.
            print('Options')
        errors = []
        for mag in ['B', 'V', 'R', 'J', 'H', 'K']:
            name = mag + 'mag'
            if np.ma.core.is_masked(catalogue[name][0]):
                errors.append(mag)
                continue
            self.set_magnitude(**{mag: catalogue[name][0]})
        if len(errors) > 0:
            warnings.warn('Magnitudes in {} were not located in NOMAD'.format(errors))
        
        
    def geocentric(self, time):
        # calculate the position of the star using parallax and proper motion
        return
    
    
    def barycentric(self, time):
        # calculate the position of the star using proper motion
        return
    
    
    def add_offset(self, da_cosdec, ddec):
        # saves an offset for the star
        
        #dadc = test_attr(da_cosdec, u.quantity.Quantity, 'd_lon_coslat')
        #dd = test_attr(ddec, u.quantity.Quantity, 'dd')
        #self.delta = SphericalCosLatDifferential(dadc, dd, 0.0*u.km)
        return
        
    
    def __str__(self):
        """String representation of the Star class
        """
        out = 'Star coordinate: RA={} +/- {}, DEC={} +/- {}\n'.format(
            self.coord.ra.to_string(u.hourangle, sep='hms', precision=5), self.errors['RA'],
            self.coord.dec.to_string(u.deg, sep='dms', precision=4), self.errors['DEC'])
        out += 'Magnitudes: '
        for mag in self.mag:
            out += '{}: {},'.format(mag, self.mag[mag])
        out += '\b\n'
        if hasattr(self, 'diameter_gaia'):
            out += 'Diameter: {}, Source: Gaia-DR2\n'.format(self.diameter_gaia)
        if hasattr(self, 'diameter_user'):
            out += 'Diameter: {}, Source: USER\n'.format(self.diameter_user)
        return out
