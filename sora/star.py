from astropy.coordinates import SkyCoord, SphericalCosLatDifferential, Distance
from astropy.coordinates import get_sun, SphericalRepresentation, SkyOffsetFrame, ICRS
from astropy.time import Time
from astropy.table import Table
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
    if 'log' in kwargs and kwargs['log']:
        print('Downloading star parameters from {}'.format(kwargs['catalog']))
    vquery = Vizier(columns=kwargs['columns'], row_limit=row_limit, timeout=600)
    if 'code' in kwargs:
        catalogue = vquery.query_constraints(catalog=kwargs['catalog'], Source=kwargs['code'])
    elif 'coord' in kwargs:
        catalogue = vquery.query_region(kwargs['coord'], radius=kwargs['radius'], catalog=kwargs['catalog'])
    else:
        raise ValueError('At least a code or coord should be given as input')
    return catalogue


def van_belle(magB=None, magV=None, magK=None):
    '''
    Determine the diameter of a star in mas using equations from van Belle (1999) 
    -- Publi. Astron. Soc. Pacific 111, 1515-1523:

    Inputs:
        magB, magV, magK: The magnitudes B, V and K of the star
    '''
    if magB is None or np.isnan(magB) or magB > 49:
        magB = np.nan
    if magV is None or np.isnan(magV) or magV > 49:
        magV = np.nan
    if magK is None or np.isnan(magK) or magK > 49:
        magK = np.nan

    def calc_diameter(a1, a2, mag):
        return 10**(a1 + a2*(mag - magK) - 0.2*mag)

    params = {'sg': {'B': [ 0.648, 0.220], 'V': [0.669, 0.223]},
              'ms': {'B': [ 0.500, 0.290], 'V': [0.500, 0.264]},
              'vs': {'B': [ 0.840, 0.211], 'V': [0.789, 0.218]}}

    mag = np.array([magB, magV])
    diameter = {}
    for st in ['sg', 'ms', 'vs']:
        diameter_s = {}
        for i,m in enumerate(['B','V']):
            diam = calc_diameter(*params[st][m], mag[i])
            if not np.isnan(diam):
                diameter_s[m] = calc_diameter(*params[st][m], mag[i])*u.mas
        if diameter_s:
            diameter[st] = diameter_s
    return diameter
    
    
def kervella(magB=None, magV=None, magK=None):
    '''
    Determine the diameter of a star in mas using equations from Kervella et. al (2004) 
    -- A&A Vol.  426, No.  1:

    Inputs:
        magB, magV, magK: The magnitudes B, V and K of the star
    '''
    if magB is None or np.isnan(magB) or magB > 49:
        magB = np.nan
    if magV is None or np.isnan(magV) or magV > 49:
        magV = np.nan
    if magK is None or np.isnan(magK) or magK > 49:
        magK = np.nan
    const1 = np.array([0.0755, 0.0535])
    const2 = np.array([0.5170, 0.5159])
    mag = np.array([magV,magB])
    vals = 10**(const1*(mag-magK)+const2-0.2*magK)
    diam = {}
    if not np.isnan(vals[0]):
        diam['V'] = vals[0]*u.mas
    if not np.isnan(vals[1]):
        diam['B'] = vals[1]*u.mas
    return diam


class Star():
    def __init__(self,**kwargs):
        '''Define a star

        Parameters:
            code (str): Gaia-DR2 Source code for searching in Vizier.
            coord (str, SkyCoord): if code is not given, coord nust have
                the coordinate to search in Vizier: 'hh mm ss.ss +dd mm ss.ss'
            nomad (bool): If true, it tries to download the magnitudes
                from NOMAD catalogue.
            log (bool): If true, it prints the downloaded information.
            local (bool): If true, it uses the given coordinate in 'coord'
                as final coordinate.
        '''
        self.__local = False
        self.mag = {}
        self.errors = {}
        self.__log = True
        if 'log' in kwargs:
            self.__log = kwargs['log']
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
            self.__searchgaia()
        nomad = kwargs.get('nomad', True)
        if nomad:
            self.__getcolors()
    
    
    def set_magnitude(self,**kwargs):
        '''
        Set the magnitudes of a star.

        Inputs:
            (band name): The magnitude for given band

        Example:
            set_magnitude(G=10)
            set_magnitude(K=15)
            set_magnitude(newband=6)
        '''
        for key in kwargs:
            mag = test_attr(kwargs[key], float, key)
            if key in self.mag:
                warnings.warn('{0} mag already defined. {0}={1} will be replaced by {0}={2}'.format(key, self.mag[key], kwargs[key]))
            self.mag[key] = kwargs[key]

            
    def set_diameter(self, diameter):
        '''
        Set a user diameter of a star in mas.

        Inputs:
            diameter = float, in mas
        '''
        self.diameter_user = diameter*u.mas

    
    def van_belle(self):
        '''
        Determine the diameter of a star in mas using equations from van Belle (1999) 
        -- Publi. Astron. Soc. Pacific 111, 1515-1523:
       '''
        return van_belle(self.mag.get('B'), self.mag.get('V'), self.mag.get('K'))
    

    def kervella(self):
        '''
        Determine the diameter of a star in mas using equations from Kervella et. al (2004) 
        -- A&A Vol.  426, No.  1:
        '''
        return kervella(self.mag.get('B'), self.mag.get('V'), self.mag.get('K'))
    
    
    def apparent_diameter(self, distance, mode='auto', log=True, **kwargs):
        """Calculate the apparent diameter of the star at given distance

        Parameters:
            distance (int, float): Object distance in AU
            mode: The way to calculate the apparent diameter
                'user': calculates using user diameter
                'gaia': calculates using Gaia diameter
                'kervella': calculates using Kervella equations
                'van_belle': calculates using van Belle equations
                'auto' (default): tries all the above methods
                    until it is able to calculate diameter.
                    The order of try is the same as shows above.
            'band': If mode is 'kervella' or 'van_belle',
                the filter must be given, 'B' or 'V'.
                In 'auto' mode, 'V' is selected.
            'star_type': If mode is 'van_belle', the star type must be given.
                It can be 'sg', 'ms' and 'vs' for 'Super Giant', 'Main
                Sequence' and 'Variable Star'. In 'auto' mode, 'sg' is selected.
            'log': If True, it prints the mode used by 'auto'.
        """
        try:
            d = distance.to(u.km)
        except:
            distance = distance*u.AU

        kwargs['band'] = kwargs.get('band', 'V')
        kwargs['star_type'] = kwargs.get('star_type', 'sg')

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
                
        if kwargs['band'] not in ['B', 'V']:
            raise KeyError('band must be informed as "B", or "V"')
            
        if mode in ['kervella', 'auto']:
            diam_kerv = self.kervella().get(kwargs['band'])
            if diam_kerv is None:
                raise ValueError('Diameter could not be calculated for given band')
            if log:
                print('Apparent diameter using Kervella et al. (2004)')
            diam = distance*np.tan(diam_kerv)
            return diam.to(u.km)
            
        if kwargs['star_type'] not in ['sg', 'ms', 'vs']:
            raise KeyError('star_type must be informed as "sg", "ms" or "vs"')
            
        if mode in ['van_belle', 'auto']:
            diam_van = self.van_belle().get(kwargs['star_type'])
            if diam_van is None:
                raise ValueError('Diameter could not be calculated using Van Belle')
            diam_van = diam_van.get(kwargs['band'])
            if diam_van is None:
                raise ValueError('Diameter could not be calculated for given band')
            if log:
                print('Apparent diameter using van Belle (1999)')
            diam = distance*np.tan(diam_van)
            return diam.to(u.km)
                        
        raise AttributeError("Star apparent diameter could not be calculated. ",
                             "Please define star diameter or B,V,K magnitudes.")
        
    
    def __searchgaia(self):
        """search for the star position in the gaia catalogue and save informations
        """
        columns = ['Source', 'RA_ICRS', 'e_RA_ICRS', 'DE_ICRS', 'e_DE_ICRS', 'Plx', 'pmRA', 'e_pmRA', 'pmDE', 'e_pmDE', 'Gmag', 'e_Gmag', 'Dup', 'Epoch', 'Rad']
        if hasattr(self, 'code'):
            catalogue = search_star(code=self.code, columns=columns, catalog='I/345/gaia2', log=self.__log)
        else:
            catalogue = search_star(coord=self.coord, columns=columns, radius=2*u.arcsec, catalog='I/345/gaia2', log=self.__log)
        if len(catalogue) == 0:
            raise ValueError('No star was found within 2 arcsec from given coordinate')
        catalogue = catalogue[0]
        if len(catalogue) > 1:
            if self.__log:
                print('{} stars were found within 2 arcsec from given coordinate.'.format(len(catalogue)))
                print('The list below is sorted by distance. Please select the correct star')
            gmag = catalogue['Gmag']
            tstars = SkyCoord(catalogue['RA_ICRS'], catalogue['DE_ICRS'])
            sep = tstars.separation(self.coord)
            k = sep.argsort()
            while True:
                t = Table()
                t['num'] = np.arange(len(tstars))+1
                t['dist(")'] = sep[k].arcsec
                t['dist(")'].format = '6.3f'
                t['Gmag'] = gmag[k].quantity.value
                t['Gmag'].format = '6.3f'
                t['RA___ICRS___DEC'] = tstars[k].to_string('hmsdms', precision=4)
                t.pprint_all()
                print('  0: Cancel')
                choice = int(input('Choose the corresponding number of the correct star: '))
                if choice in np.arange(len(k)+1):
                    break
                print('{} is not a valid choice. Please select the correct star'.format(choice))
            if choice == 0:
                raise ValueError('It was not possible to define a star')
            catalogue = catalogue[[k[choice-1]]]
        self.code = catalogue['Source'][0]
        ra = catalogue['RA_ICRS']
        dec = catalogue['DE_ICRS']
        pmra = catalogue['pmRA']
        pmde = catalogue['pmDE']
        epoch = Time(catalogue['Epoch'].quantity, format='jyear')
        try:
            distance = Distance(parallax=catalogue['Plx'].quantity, allow_negative=False)
            self.coord = SkyCoord(ra, dec, distance=distance, pm_ra_cosdec=pmra, pm_dec=pmde, obstime=epoch[0])[0]
        except:
            self.coord = SkyCoord(ra, dec, pm_ra_cosdec=pmra, pm_dec=pmde, obstime=epoch[0])[0]
        self.set_magnitude(G=catalogue['Gmag'][0])
        self.errors['RA'] = catalogue['e_RA_ICRS'][0]*u.mas
        self.errors['DEC'] = catalogue['e_DE_ICRS'][0]*u.mas
        self.errors['pmRA'] = catalogue['e_pmRA'][0]*(u.mas/u.yr)
        self.errors['pmDEC'] = catalogue['e_pmDE'][0]*(u.mas/u.yr)
        rad = catalogue['Rad'][0]
        if np.ma.core.is_masked(rad) or np.ma.core.is_masked(catalogue['Plx'][0]):
            if self.__log:
                warnings.warn('Gaia catalogue does not have star radius.')
        else:
            self.diameter_gaia = 2*np.arctan((rad*u.solRad)/distance[0]).to(u.mas)
        if self.__log:
            print('1 Gaia-DR2 star found G={}'.format(catalogue['Gmag'][0]))
            print('star coordinate at J{}: RA={} +/- {}, DEC={} +/- {}'.format(self.coord.obstime.jyear,
                  self.coord.ra.to_string(u.hourangle, sep='hms', precision=5), self.errors['RA'],
                  self.coord.dec.to_string(u.deg, sep='dms', precision=4), self.errors['DEC']))
        
        
    def __getcolors(self):
        """search for the B,V,K magnitudes of the star on Vizier and saves the result
        """
        columns = ['RAJ2000', 'DEJ2000', 'Bmag', 'Vmag', 'Rmag', 'Jmag', 'Hmag', 'Kmag']
        catalogue = search_star(coord=self.coord, columns=columns, radius=2*u.arcsec, catalog='I/297/out', log=self.__log)
        if len(catalogue) == 0:
            if self.__log:
                warnings.warn('No star was found on NOMAD that matches the star')
            return
        catalogue = catalogue[0]
        if len(catalogue) > 1:
            print('{} stars were found within 2 arcsec from given coordinate.'.format(len(catalogue)))
            print('The list below is sorted by distance. Please select the correct star')
            print('Star G mag: {}'.format(self.mag['G']))
            tstars = SkyCoord(catalogue['RAJ2000'], catalogue['DEJ2000'])
            sep = tstars.separation(self.coord)
            k = sep.argsort()
            while True:
                t = Table()
                t['num'] = np.arange(len(tstars))+1
                t['dist(")'] = sep[k].arcsec
                t['dist(")'].format = '6.3f'
                t['Bmag'] = catalogue['Bmag'][k].quantity.value
                t['Vmag'] = catalogue['Vmag'][k].quantity.value
                t['Rmag'] = catalogue['Rmag'][k].quantity.value
                t['Jmag'] = catalogue['Jmag'][k].quantity.value
                t['Hmag'] = catalogue['Hmag'][k].quantity.value
                t['Kmag'] = catalogue['Kmag'][k].quantity.value
                t['Bmag'].format = t['Vmag'].format = t['Rmag'].format = t['Jmag'].format = t['Hmag'].format = t['Kmag'].format = '6.3f'
                t['RA___ICRS___DEC'] = tstars[k].to_string('hmsdms', precision=4)
                t.pprint_all()
                print('  0: Cancel')
                choice = int(input('Choose the corresponding number of the correct star: '))
                if choice in np.arange(len(k)+1):
                    break
                print('{} is not a valid choice. Please select the correct star'.format(choice))
            if choice == 0:
                warnings.warn('No magnitudes were obtained from NOMAD')
                return
            catalogue = catalogue[[k[choice-1]]]
        errors = []
        for mag in ['B', 'V', 'R', 'J', 'H', 'K']:
            name = mag + 'mag'
            if np.ma.core.is_masked(catalogue[name][0]):
                errors.append(mag)
                continue
            self.set_magnitude(**{mag: catalogue[name][0]})
        if len(errors) > 0 and self.__log:
            warnings.warn('Magnitudes in {} were not located in NOMAD'.format(errors))
        
        
    def geocentric(self, time):
        """ calculate the position of the star using parallax and proper motion

        Parameters:
            time (float, Time):time to apply proper motion and calculate paralax.
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
            g_coord = SkyCoord(*(n_coord.cartesian.xyz + sun.cartesian.xyz), representation_type='cartesian')
            g_coord = g_coord.represent_as(SphericalRepresentation)
            g_coord = SkyCoord(g_coord.lon, g_coord.lat, g_coord.distance)

        if hasattr(self, 'offset'):
            star_frame = SkyOffsetFrame(origin=g_coord)
            new_pos = SkyCoord(lon=self.offset.d_lon_coslat, lat=self.offset.d_lat, frame=star_frame)
            return new_pos.transform_to(ICRS)
        
        return g_coord
    
    
    def barycentric(self, time):
        """ calculate the position of the star using proper motion
        
        Parameters:
            time (float, Time):time to apply proper motion.
        """
        try:
            time = Time(time)
        except:
            time = Time(time, format='jd', scale='utc')
        if np.isnan(self.coord.pm_dec):
            return self.coord
        if self.coord.distance.unit.is_unity():
            star_frame = SkyOffsetFrame(origin=self.coord)
            dt = time - self.coord.obstime
            new_pos = SkyCoord(lon=self.coord.pm_ra_cosdec*dt, lat=self.coord.pm_dec*dt, frame=star_frame)
            n_coord = new_pos.transform_to(ICRS)
        else:
            n_coord = self.coord.apply_space_motion(new_obstime=time)
        return n_coord

    def error_at(self, time):
        """Estimates star position error at time given

        Parameters:
            time (float, Time):time to apply proper motion.

        Return:
            pair of errors, in RA* and DEC
        """
        try:
            time = Time(time)
        except:
            time = Time(time, format='jd', scale='utc')
        e_ra = self.errors['RA'] + self.errors['pmRA']*np.abs(time-self.coord.obstime)
        e_dec = self.errors['DEC'] + self.errors['pmDEC']*np.abs(time-self.coord.obstime)
        return e_ra, e_dec

    def add_offset(self, da_cosdec, ddec):
        """Add an offset to star

        Parameters:
            da_cosdec (int, float):Delta_alpha_cos_delta in mas
            ddec (int, float):Delta_delta in mas
        """
        dadc = test_attr(da_cosdec, float, 'da_cosdec')
        dd = test_attr(ddec, float, 'ddec')
        self.offset = SphericalCosLatDifferential(dadc*u.mas, dd*u.mas, 0.0*u.km)

    def __str__(self):
        """String representation of the Star class
        """
        out = 'ICRS star coordinate at J{}:\nRA={} +/- {:.4f}, DEC={} +/- {:.4f}\n\n'.format(self.coord.obstime.jyear,
            self.coord.ra.to_string(u.hourangle, sep='hms', precision=5), self.errors['RA'],
            self.coord.dec.to_string(u.deg, sep='dms', precision=4), self.errors['DEC'])
        if hasattr(self, 'code'):
            out += 'Gaia-DR2 star Source ID: {}\n'.format(self.code)
        if hasattr(self, 'offset'):
            out += 'Offset Apllied: d_alpha_cos_dec = {}, d_dec = {}\n'.format(self.offset.d_lon_coslat, self.offset.d_lat)
        out += 'Magnitudes:'
        mag_out = [' {}: {:6.3f}'.format(mag, self.mag[mag]) for mag in self.mag]
        out_mag = []
        for i,mag in enumerate(mag_out):
            if i%6 == 0:
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
            out += '   ' + ','.join([' {}: {:.4f}'.format(k, v) for k,v in kerv.items()])
        vanb = self.van_belle()
        if vanb:
            out += '\nApparent diameter from van Belle (1999):'
            for key, value in vanb.items():
                out += '\n    {}:'.format(key)
                out += ','.join([' {}: {:.4f}'.format(k, v) for k,v in value.items()])
        return out
