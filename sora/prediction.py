from .star import Star
from .ephem import EphemKernel, EphemJPL, EphemPlanete
from .observer import Observer
import astropy.units as u
import astropy.constants as const
from astropy.coordinates import SkyCoord
from astropy.time import Time
from astropy.table import Table, Row
from astroquery.vizier import Vizier
import numpy as np
import warnings
import os

class Prediction(Table):
    """An Astropy Table object modified for Prediction purposes.
    """
    def __init__(self, *args, **kwargs):
        if not all([i in kwargs for i in ['time', 'coord_star', 'coord_obj','ca', 'pa', 'vel', 'dist']]):
            super().__init__(*args, names=['Epoch', 'ICRS Star Coord at Epoch', 'Geocentric Object Position', 'C/A', 'PA', 'Vel', 'Dist', 'G', 'G*', 'long', 'loct'], **kwargs)
        else:
            values = {}
            time = Time(kwargs['time'])
            size = len(time)
            values['Epoch'] = kwargs['time']; del kwargs['time']
            coord = SkyCoord(kwargs['coord_star'], unit=(u.hourangle,u.deg))
            values['ICRS Star Coord at epoch'] = kwargs['coord_star']; del kwargs['coord_star']
            values['Geocentric Object Position'] = kwargs['coord_obj']; del kwargs['coord_obj']
            values['C/A'] = kwargs['ca']; del kwargs['ca']
            values['P/A'] = kwargs['pa']; del kwargs['pa']
            values['Vel'] = kwargs['vel']; del kwargs['vel']
            values['Dist'] = kwargs['dist']; del kwargs['dist']
            if 'mag' not in kwargs.keys() and 'mag_20' not in kwargs.keys():
                raise ValueError('User must provide "mag" or "mag_20" parameters')
            if 'mag' in kwargs.keys():
                values['G'] = kwargs['mag']
                del kwargs['mag']
            else:
                values['G'] = list(['{:6.3f}'.format(float(kwargs['mag_20'][i]) - 2.5*np.log10(np.absolute(float(values['Vel'][i]))/20.0)) for i in range(size)])
            if 'mag_20' in kwargs.keys():
                values['G*'] = kwargs['mag_20']
                del kwargs['mag_20']
            else:
                values['G*'] = list(['{:6.3f}'.format(float(values['G'][i]) + 2.5*np.log10(np.absolute(float(values['Vel'][i]))/20.0)) for i in range(size)])
            if 'long' in kwargs.keys():
                values['long'] = kwargs['long']
                del kwargs['long']
            else:
                values['long'] = list(['{:3.0f}'.format(i.deg) for i in (coord.ra - time.sidereal_time('mean', 'greenwich')).wrap_at(360*u.deg)])
            if 'loct' in kwargs.keys():
                values['loct'] = kwargs['loct']
                del kwargs['loct']
            else:
                values['loct'] = list(['{}'.format((time[i] + float(values['long'][i])*(24.0/360.0)*u.hour).iso[11:16]) for i in range(size)])
            super().__init__(values, **kwargs)

    def __getitem__(self, item):
        if isinstance(item, str) and item not in self.colnames:
            col = self['Epoch']
            arr = list([i for i, c in enumerate(col) if c.startswith(item)])
            if len(arr) is not 1:
                raise KeyError('Given key is not enough to identify an unique row.')
            else:
                return self.Row(self, arr[0])
        return super().__getitem__(item)

    @classmethod
    def from_praia(cls, filename, name, **kwargs):
        from .ephem import read_obj_data
        occs = {}
        if os.path.isfile(filename) == False:
            raise IOError('File {} not found'.format(filename))
        try:
            dados = np.loadtxt(filename, skiprows=41, usecols=(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,
                                                               16, 17, 18, 19, 20, 21, 22, 25, 26, 27, 28, 29),
                dtype={'names': ('dia', 'mes', 'ano', 'hor', 'min', 'sec', 'afh', 'afm', 'afs', 'ded', 'dem', 'des', 'afh2', 'afm2', 'afs2',
                                 'ded2', 'dem2', 'des2', 'ca', 'pa', 'vel', 'delta', 'mR', 'mK', 'long', 'loct', 'ora', 'ode'),
                       'formats': ('S30', 'S30', 'S30','S30', 'S30', 'S30','S20', 'S20', 'S20','S20', 'S20', 'S20', 'S20', 'S20', 'S20',
                                   'S20', 'S20', 'S20', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'S20', 'f4', 'f4')}, ndmin=1)
        except:
            raise IOError('{} is not in PRAIA format or does not have any occultation'.format(filename))
            
        f = open(filename, 'r')
        lines = f.readlines()
        f.close()
        
        max_ca = float(lines[14].split()[-2])*u.arcsec

        ################## reading coordinates #################
        coor = np.char.array(dados['afh'], unicode=True)
        for i in ['afm', 'afs', 'ded', 'dem', 'des']:
            coor = np.core.defchararray.add(coor, ' ')
            coor = np.core.defchararray.add(coor, np.char.array(dados[i], unicode=True))
        coord = SkyCoord(coor, frame='icrs', unit=(u.hourangle, u.degree))
        coord_star = [c.to_string('hmsdms',precision=5, sep=' ') for c in coord]
        
        coor = np.char.array(dados['afh2'], unicode=True)
        for i in ['afm2', 'afs2', 'ded2', 'dem2', 'des2']:
            coor = np.core.defchararray.add(coor, ' ')
            coor = np.core.defchararray.add(coor, np.char.array(dados[i], unicode=True))
        coord = SkyCoord(coor, frame='icrs', unit=(u.hourangle, u.degree))
        coord_obj = [c.to_string('hmsdms',precision=5, sep=' ') for c in coord]

        ################### reading time ########################
        tim=np.char.array(dados['ano'], unicode=True)
        len_iso = ['-', '-', ' ', ':',':']
        arr = ['mes', 'dia', 'hor', 'min', 'sec']
        for i in np.arange(len(arr)):
            tim = np.core.defchararray.add(tim, len_iso[i]) 
            tim = np.core.defchararray.add(tim, np.char.array(dados[arr[i]], unicode=True))
        tim2 = Time(np.char.array(tim) + '000')
        time = [i.iso for i in tim2]

        ############### defining parameters #############
        ca = ['{:5.3f}'.format(i) for i in dados['ca']]
        pa = ['{:6.2f}'.format(i) for i in dados['pa']]
        vel = ['{:-6.2f}'.format(i) for i in dados['vel']]
        dist = ['{:7.3f}'.format(i) for i in dados['delta']]
        mag_20 = ['{:6.3f}'.format(i) for i in dados['mR']]
        long = ['{:3.0f}'.format(i) for i in dados['long']]
        loct = ['{:5s}'.format(i.decode()) for i in dados['loct']]
        occs['ob_off_ra'] = dados['ora']*u.mas
        occs['ob_off_de'] = dados['ode']*u.mas

        data = read_obj_data()
        radius, error_ra, error_dec = data.get(name.lower(), [0,0,0])
        radius = radius*u.km
        if 'radius' in kwargs:
            radius = kwargs['radius']*u.km
        meta = {'name': name, 'radius': radius, 'max_ca': max_ca, 'ephem': lines[17].split()[-1], 'error_ra': error_ra*1000, 'error_dec': error_dec*1000}
        return cls(time=time, coord_star=coord_star, coord_obj=coord_obj, ca=ca, pa=pa, vel=vel, mag_20=mag_20, dist=dist, long=long, loct=loct, meta=meta)

    def to_praia(self, filename):
        """ Write prediction table to PRAIA format.

        INPUT:
            filename(str): name of the file to save table
        """
        from .config import praia_occ_head
        f = open(filename, 'w')
        f.write(praia_occ_head.format(max_ca=self.meta['max_ca'].to(u.arcsec), size=len(self), ephem=self.meta.get('ephem', 'ephem')))
        for time, coord, coord_geo, ca, pa, vel, dist, mag, mag_20, long, loct in self.iterrows():
            dmag = float(mag_20)-float(mag)
            f.write("\n {} {} {}  {}  {} {}   {} {}   {}  {} {} {:5.2f} {} {:-4.1f} {:-4.1f} {:-4.1f}   {}. {}       0.0      0.0 ok g2 0    0    0    0    0".
                    format(time[8:10], time[5:7], time[:4], time[11:21].replace(':', ' '), coord[:13], coord[15:29], coord_geo[:13], coord_geo[15:29], ca, pa, vel,
                           float(dist), mag_20[:4], dmag, dmag, dmag, long, loct)
                   )
        f.close()

    def to_ow(self, ow_des, mode='append'):
        """ Write Prediction Table to OccultWatcher feeder update files.
        Tables will be saved in two files: "tableOccult_update.txt" and "LOG.dat"

        INPUT:
            ow_des (str): Occult Watcher designation for the object.
            mode (str): 'append' to append table to already existing file, default
                        'restart' to overwrite existing file.
        """
        from .config import ow_occ_head
        modes = {'append': 'a+', 'restart': 'w'}
        if mode not in modes.keys():
            raise ValueError('mode param must be "append" or "restart".')

        f = open('tableOccult_update.txt', modes[mode])
        f.write(ow_occ_head.format(name=self.meta['name'], ephem=self.meta.get('ephem', 'ephem'), max_ca=self.meta['max_ca'].to(u.arcsec),
                                   size=len(self), radius=self.meta['radius'], ow_des=ow_des))
        for time, coord, coord_geo, ca, pa, vel, dist, mag, mag_20, long, loct in self.iterrows():
            dmag = float(mag_20)-float(mag)
            f.write('{} {} {}  {}   {} {}   {} {}   {}  {} {}0 {} {} {:-4.1f}   {}. {}  {:4.0f}  {:4.0f}\n'.
                    format(time[8:10], time[5:7], time[:4], time[11:20].replace(':', ' '), coord[:13], coord[15:28], coord_geo[:13],
                           coord_geo[15:28], ca, pa, vel, dist, mag_20[:4], dmag, long, loct, self.meta.get('error_ra', 0),
                           self.meta.get('error_dec', 0)))
        f.write(' '+'-'*148+'\n')
        f.close()

        f = open('LOG2.dat', modes[mode])
        t = Time.now()
        for time in Time(self['Epoch']):
            t0 = Time(time.isot.split('T')[0] + ' 00:00:00.0')
            dt = (time-t0).jd*24
            f.write('{} {:5s} {}-{:06.3f}\n'.format(t.isot[:-7], ow_des, time.isot.split('T')[0], dt))
        f.close()


def occ_params(star, ephem, time):
    """ Calculates the parameters of the occultation, as instant, CA, PA.
        
    Parameters:
    star (Star): The coordinate of the star in the same frame as the ephemeris.
    It must be a Star object.
    ephem (Ephem): Ephemeris. It must be an Ephemeris object.
    
    Return:
    instant of CA (Time): Instant of Closest Approach
    CA (arcsec): Distance of Closest Approach
    PA (deg): Position Angle at Closest Approach
    """
    
    delta_t = 0.05
    
    if type(star) != Star:
        raise ValueError('star must be a Star object')
    if type(ephem) not in [EphemKernel, EphemJPL, EphemPlanete]:
        raise ValueError('ephem must be an Ephemeris object')
        
    tt = time + np.arange(-600, 600, delta_t)*u.s
    coord = star.geocentric(tt[0])
    if type(ephem) == EphemPlanete:
        ephem.fit_d2_ksi_eta(coord, log=False)
    
    if type(ephem) == EphemJPL:
        tt = time + np.arange(-600, 600, 4)*u.s
        ksi, eta = ephem.get_ksi_eta(tt, coord)
        dd = np.sqrt(ksi*ksi+eta*eta)
        min = np.argmin(dd)
        tt = tt[min] + np.arange(-8, 8, 0.05)*u.s
    
    ksi, eta = ephem.get_ksi_eta(tt, coord)
    dd = np.sqrt(ksi*ksi+eta*eta)
    min = np.argmin(dd)
    
    if type(ephem) == EphemPlanete:
        dist = ephem.ephem[int(len(ephem.time)/2)].distance
    else:
        dist = ephem.get_position(time).distance
    
    ca = np.arcsin(dd[min]*u.km/dist).to(u.arcsec)
    
    pa = (np.arctan2(-ksi[min],-eta[min])*u.rad).to(u.deg)
    if pa < 0*u.deg:
        pa = pa + 360*u.deg
    
    dksi = ksi[min+1]-ksi[min]
    deta = eta[min+1]-eta[min]
    vel = np.sqrt(dksi**2 + deta**2)/delta_t
    vel = -vel*np.sign(dksi)*(u.km/u.s)
    
    return tt[min], ca, pa, vel, dist.to(u.AU)


def prediction(ephem, time_beg, time_end, mag_lim=None, interv=60, divs=1, sigma=1):
    """ Predicts stellar occultations

    Parameters:
    ephem (Ephem): Ephemeris. It must be an Ephemeris object.
    time_beg (Time): Initial time for prediction
    time_beg (Time): Final time for prediction
    mag_lim (int,float): Faintest Gmag for search
    interv (int, float): interval, in seconds, of ephem times for search
    divs (int,float): interal, in deg, for max search of stars

    Return:
    occ_params (Table): Astropy Table with the occultation params for each event
    """
    # generate ephemeris
    if type(ephem) != EphemKernel:
        raise TypeError('At the moment prediction only works with EphemKernel')
    print("Generating Ephemeris ...")
    time_beg = Time(time_beg)
    time_end = Time(time_end)
    dt = np.arange(0, (time_end-time_beg).sec, interv)*u.s
    t = time_beg + dt
    coord = ephem.get_position(t)

    # define catalogue parameters
    kwds = {}
    kwds['columns'] = ['Source', 'RA_ICRS', 'DE_ICRS']
    kwds['row_limit']=10000000
    kwds['timeout']=600
    if mag_lim:
        kwds['column_filters']={"Gmag":"<{}".format(mag_lim)}
    vquery = Vizier(**kwds)

    # determine suitable divisions for star search
    radius = ephem.radius + const.R_earth
    mindist = np.arcsin(radius/coord[0].distance)
    divisions = []
    n=0
    while True:
        dif = coord.separation(coord[n])
        k = np.where(dif < divs*u.deg)[0]
        l = np.where(k[1:]-k[:-1] > 1)[0]
        l = np.append(l,len(k)-1)
        m = np.where(k[l] - n > 0)[0].min()
        k = k[l][m]
        divisions.append([n,k])
        if k == len(coord)-1:
            break
        n = k

    print('Ephemeris was split in {} parts for better search of stars'.format(len(divisions)))

    # makes predictions for each division
    occs = []
    for i,vals in enumerate(divisions):
        print('\nSearching occultations in part {}/{}'.format(i+1,len(divisions)))
        nt = t[vals[0]:vals[1]]
        ncoord = coord[vals[0]:vals[1]]
        ra = np.mean([ncoord.ra.min().deg,ncoord.ra.max().deg])
        dec = np.mean([ncoord.dec.min().deg,ncoord.dec.max().deg])
        width = ncoord.ra.max() - ncoord.ra.min() + 2*mindist
        height = ncoord.dec.max() - ncoord.dec.min() + 2*mindist
        pos_search = SkyCoord(ra*u.deg, dec*u.deg)

        print('Downloading stars ...')
        catalogue = vquery.query_region(pos_search, width=width, height=height, catalog='I/345/gaia2')
        print('Identifying occultations ...')
        if len(catalogue) == 0:
            continue
        catalogue = catalogue[0]
        stars = SkyCoord(catalogue['RA_ICRS'], catalogue['DE_ICRS'])
        idx, d2d, d3d = stars.match_to_catalog_sky(ncoord)

        dist = np.arcsin(radius/ncoord[idx].distance) + sigma*np.max([ephem.error_ra.value,ephem.error_dec.value])*u.arcsec
        k = np.where(d2d < dist)[0]
        for ev in k:
            star = Star(code=catalogue['Source'][ev], log=False)
            pars = [star.code, star.geocentric(nt[idx][ev]), star.mag['G']]
            pars = np.hstack((pars, occ_params(star, ephem, nt[idx][ev])))
            occs.append(pars)

    if not occs:
        warnings.warn('No stellar occultation was found')
        return Prediction(meta=meta)
    # create astropy table with the params
    meta = {'name': ephem.name, 'time_beg': time_beg, 'time_end': time_end, 'maglim': mag_lim, 'max_ca': dist.max(),
            'radius':ephem.radius.to(u.km).value, 'error_ra': ephem.error_ra.to(u.mas).value, 'error_dec': ephem.error_dec.to(u.mas).value}
    occs2 = np.transpose(occs)
    time = Time(occs2[3])
    geocentric = ephem.get_position(time)
    k = np.argsort(time)
    source = occs2[0][k]
    coord = [i.to_string('hmsdms',precision=5, sep=' ') for i in occs2[1][k]]
    coord_geo = [i.to_string('hmsdms',precision=5, sep=' ') for i in geocentric[k]]
    mags = ['{:6.3f}'.format(i) for i in occs2[2][k]]
    mags_20 = ['{:6.3f}'.format(occs2[2][i] + 2.5*np.log10(np.absolute(occs2[6][i].value)/20.0)) for i in k]
    time = [i.iso for i in time[k]]
    ca = ['{:5.3f}'.format(i.value) for i in occs2[4][k]]
    pa = ['{:6.2f}'.format(i.value) for i in occs2[5][k]]
    vel = ['{:-6.2f}'.format(i.value) for i in occs2[6][k]]
    dist = ['{:7.3f}'.format(i.value) for i in occs2[7][k]]
    t = Prediction(time=time, coord_star=coord, coord_obj=coord_geo, ca=ca, pa=pa, vel=vel, mag=mags, mag_20=mags_20, dist=dist, meta=meta)
    return t
