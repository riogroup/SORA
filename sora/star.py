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
    print('Downloading star parameters from Gaia-DR2')
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

class Star():
    '''
    Docstring
    Define a star
    '''
    def __init__(self,**kwargs):
        self.__local = False
        self.mags = {}
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
            if key in self.mag:
                warnings.warn('{0} mag already defined. {1} will be replaced by {2}'.format(key, self.mag[key], kwargs[key]))
            self.mag[key] = kwargs[key]

    def set_diameter(self,star_diameter):
        '''
        Set the diameter of a star in mas.
        Inputs:
        star_diameter = float, in mas
        '''
        self.diameter = star_diameter

    def calc_diameter(self,obs_filter='V',star_type='sg'):
        '''
        Determine the diameter of a star in mas using equations from van Belle (1999) 
        -- Publi. Astron. Soc. Pacific 111, 1515-1523:
        Inputs:
        obs_filter = str; should be Visual (V) or Blue (B)
        star_type = str; should be Super-giant (gs), Main-sequence (ms) or Variable (vs)
        '''
        if (self.magV == None):
            print('First you need to set the magnitudes B, V, K')
            return
        if (obs_filter not in ['V','B']):
            print('obs_filter should be Visual (V) or Blue (B)')
            return
        if (star_type not in ['sg','ms','vs']):
            print('star_type should be Super-giant (gs), Main-sequence (ms) or Variable (vs)')
            return
        return van_belle(self.magB, self.magV, self.magK)[star_type][obs_filter]
    
    def apparent_radius(self, distance):
        # calculate the apparent radius of the star at given distance
        return
    
    def __searchgaia2(self):
        """search for the star position in the gaia catalogue and save informations
        """
        columns = ['Source', 'RA_ICRS', 'e_RA_ICRS', 'DE_ICRS', 'e_DE_ICRS', 'Plx', 'pmRA', 'e_pmRA', 'pmDE', 'e_pmDE', 'Gmag', 'e_Gmag', 'Dup', 'Epoch', 'Rad']
        if hasattr(self, 'code'):
            catalogue = search_star(code=self.code, columns=columns, catalog='I/345/gaia2')[0]
        else:
            catalogue = search_star(coord=self.coord, columns=columns, radius=2*u.arcsec, catalog='I/345/gaia2')[0]
        if len(catalogue) == 1:
            print('1 star found')
            self.code = catalogue['Source']
            ra = catalogue['RA_ICRS']
            dec = catalogue['DE_ICRS']
            distance = Distance(parallax=catalogue['Plx'].quantity, allow_negative=True)
            pmra = catalogue['pmRA']
            pmde = catalogue['pmDE']
            epoch = Time(catalogue['Epoch'].quantity, format='jyear')
            self.coord = SkyCoord(ra, dec, distance=distance, pm_ra_cosdec=pmra, pm_dec=pmde, obstime=epoch)[0]
            self.set_magnitude(G=catalogue['Gmag'][0])
            self.errors['RA'] = catalogue['e_RA_ICRS'][0]
            self.errors['DEC'] = catalogue['e_DE_ICRS'][0]
            self.errors['pmRA'] = catalogue['e_pmRA'][0]
            self.errors['pmDEC'] = catalogue['e_pmDEC'][0]
            rad = catalogue['Rad'][0]
            if np.ma.core.is_masked(rad):
                warnings.warn('Gaia star does not have Radius, please define [B, V, K] magnitudes.')
            else:
                self.radius = rad*u.solRad
        else:
            ## pegar todas as estrelas e colocar como opcao pro usuario escolher.
            warnings.warn('{} stars found in the region searched'.format(len(catalogue)))
        
    def __getcolors(self):
        # search for the B,V,K magnitudes of the star on Vizier and saves the result
        columns = ['RAJ2000', 'DEJ2000', 'Bmag', 'Vmag', 'Rmag', 'Jmag', 'Hmag', 'Kmag']
        catalogue = search_star(coord=self.coord, columns=columns, radius=2*u.arcsec, catalog='	I/297/out')[0]
        if len(catalogue) == 0:
            raise Error('No star was found on NOMAD that matches the star')
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
        if len(errors) > 1:
            warnings.warn('Magnitudes in {} were not located in NOMAD'.format())
        
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
        # return what it is to be printed
        return ''
