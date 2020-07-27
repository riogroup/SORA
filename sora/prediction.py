from .star import Star
from .ephem import EphemKernel, EphemJPL, EphemPlanete
from sora.config import input_tests
import astropy.units as u
import astropy.constants as const
from astropy.coordinates import SkyCoord, EarthLocation, Angle, get_sun
from astropy.coordinates import get_moon, GCRS, ITRS, SkyOffsetFrame
from astropy.time import Time
from astropy.table import Table, Row, Column
from astroquery.vizier import Vizier
import numpy as np
import warnings
import os
import matplotlib.pyplot as plt
import glob


class PredictRow(Row):
    """ An Astropy Row object modified for Prediction purposes.
    """
    def plot_occ_map(self, **kwargs):
        """
        Parameters:
            radius (int,float): The radius of the shadow. If not given it uses saved value

            nameimg (str): Change the name of the image saved.
            path (str): Path to a directory where to save map.
            resolution (int): Cartopy feature resolution. "1" means a resolution of "10m",
                "2" a resolution of "50m" and "3" a resolution of "100m". Default = 2
            states (bool): True to plot the state division of the countries. The states of
                some countries will only be shown depending on the resolution.
            zoom (int, float): Zooms in or out of the map.
            centermap_geo (list): Center the map given coordinates in longitude and latitude.
                It must be a list with two numbers. Default=None.
            centermap_delta (list): Displace the center of the map given displacement
                in X and Y, in km. It must be a list with two numbers. Default=None.
            centerproj (list): Rotates the Earth to show occultation with the center
                projected at a given longitude and latitude. It must be a list with two numbers
            labels (bool): Plots text above and below the map with the occultation parameters.
                Default=True.
            meridians (int): Plots lines representing the meridians for given interval. Default=30 deg
            parallels (int): Plots lines representing the parallels for given interval. Default=30 deg
            sites (dict): Plots site positions in map. It must be a python dictionary where the key is
                the name of the site, and the value is a list with longitude, latitude, delta_x,
                delta_y and color. delta_x and delta_y are displacement, in km, from the point
                of the site in the map and the name. color is the color of the point.
            site_name (bool): If True, it prints the name of the sites given, else it plots only the points
            countries (dict): Plots the names of countries. It must be a python dictionary where the key
                is the name of the country and the value is a list with longitude and latitude
                of the lower left part of the text.
            offset (list): applies an offset to the ephemeris, calculating new CA and instant of CA.
                It is a pair of delta_RA*cosDEC and delta_DEC.
            mapstyle (int): Define the color style of the map. 1 is the default black and white scale.
                2 is a colored map.
            error (int,float): Ephemeris error in mas. It plots a dashed line representing radius + error.
            ercolor (str): Changes the color of the lines of the error bar.
            ring (int,float): It plots a dashed line representing the location of a ring.
                It is given in km, from the center.
            rncolor (str): Changes the color of ring lines.
            atm (int,float): plots a dashed line representing the location of an atmosphere.
                It is given in km, from the center.
            atcolor (str): Changes the color of atm lines.
            chord_delta (list): list with distances from center to plot chords
            chord_geo (2d-list): list with pairs of coordinates to plot chords
            chcolor (str): color of the line of the chords. Default: grey
            heights (list): plots a circular dashed line showing the locations where the observer
                would observe the occultation at a given height above the horizons.
                This must be a list.
            hcolor (str): Changes the color of the height lines.
            mapsize (list): The size of figure, in cm. It must be a list with two values.
                Default = [46.0, 38.0].
            cpoints (int,float): Interval for the small points marking the center of shadow,
                in seconds. Default=60.
            ptcolor (str): Change the color of the center points.
            alpha (float): The transparency of the night shade, where 0.0 is full transparency
                and 1.0 is full black. Default = 0.2.
            fmt (str): The format to save the image. It is parsed directly by matplotlib.pyplot.
                Default = 'png'
            dpi (int): "Dots per inch". It defines the quality of the image. Default = 100.
            lncolor (str): Changes the color of the line that represents the limits of the shadow over Earth.
            outcolor (str): Changes the color of the lines that represents the limits of the shadow outside Earth
            nscale (int,float): Arbitrary scale for the size of the name of the site.
            cscale (int,float): Arbitrary scale for the name of the country.
            sscale (int,float): Arbitrary scale for the size of point of the site.
            pscale (int,float): Arbitrary scale for the size of the points that represent the center of the shadow
            arrow (bool): If true, it plots the arrow with the occultation direction.

            Comment: Only one of centermap_geo and centermap_delta can be given
        """
        radius = kwargs.pop('radius', self.meta['radius'])
        plot_occ_map(self.meta['name'], radius, coord=self['ICRS Star Coord at Epoch'], time=self['Epoch'],
                     ca=float(self['C/A']), pa=float(self['P/A']), vel=float(self['Vel']), dist=float(self['Dist']),
                     mag=float(self['G*']), longi=float(self['long']), **kwargs)


class PredictionTable(Table):
    """ An Astropy Table object modified for Prediction purposes.
    """
    Row = PredictRow

    def __init__(self, *args, **kwargs):
        if 'time' in kwargs:
            values = {}
            time = Time(kwargs['time'])
            time.delta_ut1_utc = 0.0
            time.format = 'iso'
            if time.isscalar:
                time = Time([time])
            values['Epoch'] = Column(time)

            def coord_fmt(value):
                return value.to_string('hmsdms', precision=5, sep=' ')

            coord = SkyCoord(kwargs['coord_star'], unit=(u.hourangle, u.deg))
            values['ICRS Star Coord at Epoch'] = Column(coord, format=coord_fmt)
            try:
                coord_geo = SkyCoord(kwargs['coord_obj'])
            except:
                coord_geo = SkyCoord(kwargs['coord_obj'], unit=(u.hourangle, u.deg))
            values['Geocentric Object Position'] = Column(coord_geo, format=coord_fmt)
            values['C/A'] = Column(kwargs['ca'], format='5.3f', unit='arcsec')
            values['P/A'] = Column(kwargs['pa'], format='6.2f', unit='deg')
            values['Vel'] = Column(kwargs['vel'], format='-6.2f', unit='km/s')
            values['Dist'] = Column(kwargs['dist'], format='7.3f', unit='AU')
            del(kwargs['time'], kwargs['coord_star'], kwargs['coord_obj'], kwargs['ca'],
                kwargs['pa'], kwargs['vel'], kwargs['dist'])
            if 'mag' not in kwargs.keys() and 'mag_20' not in kwargs.keys():
                raise ValueError('User must provide "mag" or "mag_20" parameters')
            if 'mag' in kwargs.keys():
                values['G'] = Column(kwargs['mag'], format='6.3f', unit='mag')
                del kwargs['mag']
            else:
                values['G'] = kwargs['mag_20'] - 2.5*np.log10(np.absolute(values['Vel'])/20.0)
                values['G'].unit = 'mag'
                values['G'].format = '6.3f'
            if 'mag_20' in kwargs.keys():
                values['G*'] = Column(kwargs['mag_20'], format='6.3f', unit='mag')
                del kwargs['mag_20']
            else:
                values['G*'] = values['G'] + 2.5*np.log10(np.absolute(values['Vel'])/20.0)
                values['G*'].unit = 'mag'
                values['G*'].format = '6.3f'
            if 'long' in kwargs.keys():
                values['long'] = Column(kwargs['long'], unit='deg', format='3.0f')
                del kwargs['long']
            else:
                values['long'] = Column((coord.ra - time.sidereal_time('mean', 'greenwich')).wrap_at(360*u.deg),
                                        unit='deg', format='3.0f')
            if 'loct' in kwargs.keys():
                values['loct'] = Column(kwargs['loct'], unit='hh:mm')
                del kwargs['loct']
            else:
                longi = (coord.ra - time.sidereal_time('mean', 'greenwich')).wrap_at(360*u.deg)
                ntime = time + longi.hour*u.hour
                values['loct'] = Column(['{}'.format(t.iso[11:16]) for t in ntime], unit='hh:mm')
            moon_pos = get_moon(time)
            values['M-G-T'] = Column(moon_pos.separation(coord), unit='deg', format='3.0f')
            sun_pos = get_sun(time)
            values['S-G-T'] = Column(sun_pos.separation(coord), unit='deg', format='3.0f')
            if 'source' in kwargs.keys():
                values['GAIA-DR2 Source ID'] = Column(kwargs['source'])
                del kwargs['source']
            else:
                values['GAIA-DR2 Source ID'] = Column(np.repeat('', len(time)))
            super().__init__(values, **kwargs)
        else:
            super().__init__(*args, **kwargs)

    def __itens_by_epoch(self, date):
        """ Gets item list for all occultations that matches the given date

        Parameters:
            date (str): date to match

        Returns:
            item (list): the list of occultations that matches the date
        """
        col = self['Epoch']
        arr = list([i for i, c in enumerate(col) if c.iso.startswith(date)])
        return arr

    def __getitem__(self, item):
        """ The redefinition of __getitem__ allows for selecting prediction based on the ISO date of the event
        """
        if isinstance(item, str) and item not in self.colnames:
            arr = self.__itens_by_epoch(item)
            if len(arr) == 0:
                raise KeyError('No prediction corresponds to time "{}"'.format(item))
            elif len(arr) > 1:
                return self[arr]
            else:
                return self.Row(self, arr[0])
        return super().__getitem__(item)

    @classmethod
    def from_praia(cls, filename, name, **kwargs):
        """ Creates a PredictionTable Table reading from a PRAIA table

        Parameters:
            filename (str): path to the PRAIA table file.
            name (str): Name of the Object of the prediction.
            radius (int,float): Object radius, in km. (not required)
                If not given it's searched in online database.
                If not found online, the default is set to zero.

        Returns:
            A PredictionTable
        """
        from .ephem import read_obj_data
        if not os.path.isfile(filename):
            raise IOError('File {} not found'.format(filename))
        input_tests.check_kwargs(kwargs, allowed_kwargs=['radius'])
        try:
            dados = np.loadtxt(
                filename, skiprows=41,
                usecols=(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 25, 26, 27, 28, 29),
                dtype={'names': ('dia', 'mes', 'ano', 'hor', 'min', 'sec', 'afh', 'afm', 'afs', 'ded',
                                 'dem', 'des', 'afh2', 'afm2', 'afs2', 'ded2', 'dem2', 'des2', 'ca',
                                 'pa', 'vel', 'delta', 'mR', 'mK', 'long', 'loct', 'ora', 'ode'),
                       'formats': ('S30', 'S30', 'S30', 'S30', 'S30', 'S30', 'S20', 'S20', 'S20',
                                   'S20', 'S20', 'S20', 'S20', 'S20', 'S20', 'S20', 'S20', 'S20',
                                   'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'S20', 'f4', 'f4')}, ndmin=1)
        except:
            raise IOError('{} is not in PRAIA format or does not have any occultation'.format(filename))

        f = open(filename, 'r')
        lines = f.readlines()
        f.close()

        max_ca = float(lines[14].split()[-2])*u.arcsec

        # reading coordinates
        coor = np.char.array(dados['afh'], unicode=True)
        for i in ['afm', 'afs', 'ded', 'dem', 'des']:
            coor = np.core.defchararray.add(coor, ' ')
            coor = np.core.defchararray.add(coor, np.char.array(dados[i], unicode=True))
        coord_star = SkyCoord(coor, frame='icrs', unit=(u.hourangle, u.degree))

        coor = np.char.array(dados['afh2'], unicode=True)
        for i in ['afm2', 'afs2', 'ded2', 'dem2', 'des2']:
            coor = np.core.defchararray.add(coor, ' ')
            coor = np.core.defchararray.add(coor, np.char.array(dados[i], unicode=True))
        coord_obj = SkyCoord(coor, frame='icrs', unit=(u.hourangle, u.degree))

        # reading time
        tim = np.char.array(dados['ano'], unicode=True)
        len_iso = ['-', '-', ' ', ':', ':']
        arr = ['mes', 'dia', 'hor', 'min', 'sec']
        for i in np.arange(len(arr)):
            tim = np.core.defchararray.add(tim, len_iso[i])
            tim = np.core.defchararray.add(tim, np.char.array(dados[arr[i]], unicode=True))
        time = Time(np.char.array(tim) + '000')

        # defining parameters
        data = read_obj_data()
        radius, error_ra, error_dec = data.get(name.lower(), [0, 0, 0])
        radius = kwargs.get('radius', radius)*u.km
        meta = {'name': name, 'radius': radius, 'max_ca': max_ca, 'ephem': lines[17].split()[-1],
                'error_ra': error_ra*1000, 'error_dec': error_dec*1000}
        return cls(time=time, coord_star=coord_star, coord_obj=coord_obj, ca=dados['ca'], pa=dados['pa'],
                   vel=dados['vel'], mag_20=dados['mR'], dist=dados['delta'],
                   long=dados['long'], loct=dados['loct'], meta=meta)

    def to_praia(self, filename):
        """ Writes PredictionTable to PRAIA format.

        Parameters:
            filename(str): name of the file to save table
        """
        from sora.config.variables import praia_occ_head
        f = open(filename, 'w')
        f.write(praia_occ_head.format(max_ca=self.meta['max_ca'].to(u.arcsec), size=len(self),
                                      ephem=self.meta.get('ephem', 'ephem')))
        for time, coord, coord_geo, ca, pa, vel, dist, mag, mag_20, longi, loct, md, sd, source in self.iterrows():
            dmag = mag_20-mag
            f.write("\n {} {} {}  {}  {}   {}   {:5.3f}  {:6.2f} {:-6.2f} {:5.2f} {:4.1f} "
                    "{:-4.1f} {:-4.1f} {:-4.1f}   {:3.0f}. {}       0.0      0.0 ok g2 0    0    0    0    0".
                    format(time.iso[8:10], time.iso[5:7], time.iso[:4], time.iso[11:21].replace(':', ' '),
                           coord.to_string('hmsdms', precision=4, sep=' '), coord_geo.to_string('hmsdms', precision=4, sep=' '),
                           ca, pa, vel, dist, mag_20, dmag, dmag, dmag, longi, loct))
        f.close()

    def to_ow(self, ow_des, mode='append'):
        """ Writes PredictionTable to OccultWatcher feeder update file format.
            Tables will be saved in two files: "tableOccult_update.txt" and "LOG.dat"

        Parameters:
            ow_des (str): Occult Watcher designation for the object.
            mode (str): 'append' to append table to already existing file, default
                'restart' to overwrite existing file.
        """
        from sora.config.variables import ow_occ_head
        modes = {'append': 'a+', 'restart': 'w'}
        if mode not in modes.keys():
            raise ValueError('mode param must be "append" or "restart".')

        f = open('tableOccult_update.txt', modes[mode])
        if self.meta['radius'] == 0.0:
            warnings.warn('radius is 0.0, please update the value manually in "tableOccult_update.txt" '
                          'before submitting the file')
        if (self.meta['error_ra'] == 0.0) or (self.meta['error_ra'] == 0.0):
            warnings.warn('error_ra and/or error_dec is 0.0, please update the value manually in '
                          '"tableOccult_update.txt" before submitting the file')
        f.write(ow_occ_head.format(name=self.meta['name'], ephem=self.meta.get('ephem', 'ephem'),
                                   max_ca=self.meta['max_ca'].to(u.arcsec), size=len(self),
                                   radius=self.meta['radius'], ow_des=ow_des))
        for time, coord, coord_geo, ca, pa, vel, dist, mag, mag_20, longi, loct, md, sd, source in self.iterrows():
            dmag = mag_20-mag
            f.write('{} {} {}  {}   {}   {}   {:5.3f}  {:6.2f} {:-7.3f} {:7.3f} '
                    '{:4.1f} {:-4.1f}   {:3.0f}. {}  {:4.0f}  {:4.0f}\n'.
                    format(time.iso[8:10], time.iso[5:7], time.iso[:4], time.iso[11:20].replace(':', ' '),
                           coord.to_string('hmsdms', precision=4, sep=' ')[:-1],
                           coord_geo.to_string('hmsdms', precision=4, sep=' ')[:-1],
                           ca, pa, vel, dist, mag_20, dmag, longi, loct,
                           self.meta.get('error_ra', 0), self.meta.get('error_dec', 0))
                    )
        f.write(' '+'-'*148+'\n')
        f.close()

        f = open('LOG.dat', modes[mode])
        t = Time.now()
        for time in self['Epoch']:
            t0 = Time(time.isot.split('T')[0] + ' 00:00:00.0')
            dt = (time-t0).jd*24
            f.write('{} {:5s} {}-{:06.3f}\n'.format(t.isot[:-7], ow_des, time.isot.split('T')[0], dt))
        f.close()

    def plot_occ_map(self, **kwargs):
        """
        Parameters:
            radius (int,float): The radius of the shadow. If not given it uses saved value

            nameimg (str): Change the name of the image saved.
            path (str): Path to a directory where to save map.
            resolution (int): Cartopy feature resolution. "1" means a resolution of "10m",
                "2" a resolution of "50m" and "3" a resolution of "100m". Default = 2
            states (bool): True to plot the state division of the countries. The states of
                some countries will only be shown depending on the resolution.
            zoom (int, float): Zooms in or out of the map.
            centermap_geo (list): Center the map given coordinates in longitude and latitude.
                It must be a list with two numbers. Default=None.
            centermap_delta (list): Displace the center of the map given displacement
                in X and Y, in km. It must be a list with two numbers. Default=None.
            centerproj (list): Rotates the Earth to show occultation with the center
                projected at a given longitude and latitude. It must be a list with two numbers
            labels (bool): Plots text above and below the map with the occultation parameters.
                Default=True.
            meridians (int): Plots lines representing the meridians for given interval. Default=30 deg
            parallels (int): Plots lines representing the parallels for given interval. Default=30 deg
            sites (dict): Plots site positions in map. It must be a python dictionary where the key is
                the name of the site, and the value is a list with longitude, latitude, delta_x,
                delta_y and color. delta_x and delta_y are displacement, in km, from the point
                of the site in the map and the name. color is the color of the point.
            site_name (bool): If True, it prints the name of the sites given, else it plots only the points
            countries (dict): Plots the names of countries. It must be a python dictionary where the key
                is the name of the country and the value is a list with longitude and latitude
                of the lower left part of the text.
            offset (list): applies an offset to the ephemeris, calculating new CA and instant of CA.
                It is a pair of delta_RA*cosDEC and delta_DEC.
            mapstyle (int): Define the color style of the map. 1 is the default black and white scale.
                2 is a colored map.
            error (int,float): Ephemeris error in mas. It plots a dashed line representing radius + error.
            ercolor (str): Changes the color of the lines of the error bar.
            ring (int,float): plots a dashed line representing the location of a ring.
                It is given in km, from the center.
            rncolor (str): Changes the color of ring lines.
            atm (int,float): plots a dashed line representing the location of an atmosphere.
                It is given in km, from the center.
            atcolor (str): Changes the color of atm lines.
            chord_delta (list): list with distances from center to plot chords
            chord_geo (2d-list): list with pairs of coordinates to plot chords
            chcolor (str): color of the line of the chords. Default: grey
            heights (list): It plots a circular dashed line showing the locations where the observer
                would observe the occultation at a given height above the horizons.
                This must be a list.
            hcolor (str): Changes the color of the height lines.
            mapsize (list): The size of figure, in cm. It must be a list with two values.
                Default = [46.0, 38.0].
            cpoints (int,float): Interval for the small points marking the center of shadow,
                in seconds. Default=60.
            ptcolor (str): Change the color of the center points.
            alpha (float): The transparency of the night shade, where 0.0 is full transparency
                and 1.0 is full black. Default = 0.2.
            fmt (str): The format to save the image. It is parsed directly by matplotlib.pyplot.
                Default = 'png'
            dpi (int): "Dots per inch". It defines the quality of the image. Default = 100.
            lncolor (str): Changes the color of the line that represents the limits of the shadow over Earth.
            outcolor (str): Changes the color of the lines that represents the limits of the shadow outside Earth
            nscale (int,float): Arbitrary scale for the size of the name of the site.
            cscale (int,float): Arbitrary scale for the name of the country.
            sscale (int,float): Arbitrary scale for the size of point of the site.
            pscale (int,float): Arbitrary scale for the size of the points that represent the center of the shadow
            arrow (bool): If true, it plots the arrow with the occultation direction.

            Comment: Only one of centermap_geo and centermap_delta can be given
        """
        for i in range(len(self)):
            self[i].plot_occ_map(**kwargs)

    def remove_occ(self, date):
        """ Removes stellar occultations from table

        Parameters:
            date (str,list): Date or list of dates of the occultation to be removed.
                The dates mut be as shown in the 'Epoch' column. If the date is not
                complete, the function will select all occultations that matches the given string.
                For instance, date='2020-06' will remove all occultations from the month of June 2020.
        """
        if type(date) == str:
            date = [date]
        itens = list(set([y for d in date for y in self.__itens_by_epoch(d)]))
        self.remove_rows(itens)

    def keep_from_selected_images(self, path='.'):
        """ Keeps predictions which images were not deleted in given path
            This function uses the name of the images to identify predictions.
            The name must be the automatic one generated by plot_occ_map().
            The format of the image is not relevant.

        Parameters:
            path (str): path where images are located
        """
        itens = []
        for i, tca in enumerate(self['Epoch']):
            name = '{}_{}'.format(self.meta['name'], tca.isot)
            if not glob.glob(os.path.join(path, name)+'*'):
                itens.append(i)
        self.remove_rows(itens)


def occ_params(star, ephem, time):
    """ Calculates the parameters of the occultation, as instant, CA, PA.

    Parameters:
        star (Star): The coordinate of the star in the same reference frame as the ephemeris.
            It must be a Star object.
        ephem (Ephem): object ephemeris. It must be an Ephemeris object.

    Returns:
        instant of CA (Time): Instant of Closest Approach
        CA (arcsec): Distance of Closest Approach
        PA (deg): Position Angle at Closest Approach
        vel (km/s): Velocity of the occultation
        dist (AU): the object geocentric distance.
    """

    delta_t = 0.05

    if type(star) != Star:
        raise ValueError('star must be a Star object')
    if type(ephem) not in [EphemKernel, EphemJPL, EphemPlanete]:
        raise ValueError('ephem must be an Ephemeris object')

    time = Time(time)
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

    pa = (np.arctan2(ksi[min], eta[min])*u.rad).to(u.deg)
    if pa < 0*u.deg:
        pa = pa + 360*u.deg

    dksi = ksi[min+1]-ksi[min]
    deta = eta[min+1]-eta[min]
    vel = np.sqrt(dksi**2 + deta**2)/delta_t
    vel = vel*np.sign(dksi)*(u.km/u.s)

    return tt[min], ca, pa, vel, dist.to(u.AU)


def prediction(ephem, time_beg, time_end, mag_lim=None, step=60, divs=1, sigma=1, log=True):
    """ Predicts stellar occultations

    Parameters:
        ephem (Ephem): object ephemeris. It must be an Ephemeris object.
        time_beg (str,Time): Initial time for prediction
        time_beg (str,Time): Final time for prediction
        mag_lim (int,float): Faintest Gmag for search
        step (int, float): step, in seconds, of ephem times for search
        divs (int): number of regions the ephemeris will be splitted for better search of occultations
        sigma (int,float): ephemeris error sigma for search off-Earth.
        log (bool): To show what is being done at the moment.

    Returns:
        predict (PredictionTable): PredictionTable with the occultation params for each event
    """
    # generate ephemeris
    if type(ephem) is not EphemKernel:
        raise TypeError('At the moment prediction only works with EphemKernel')
    time_beg = Time(time_beg)
    time_end = Time(time_end)
    dt = np.arange(0, (time_end-time_beg).sec, step)*u.s
    t = time_beg + dt

    # define catalogue parameters
    kwds = {}
    kwds['columns'] = ['Source', 'RA_ICRS', 'DE_ICRS']
    kwds['row_limit'] = 10000000
    kwds['timeout'] = 600
    if mag_lim:
        kwds['column_filters'] = {"Gmag": "<{}".format(mag_lim)}
    vquery = Vizier(**kwds)

    # determine suitable divisions for star search
    radius = ephem.radius + const.R_earth

    divisions = np.array_split(np.arange(len(t)), divs)

    if log:
        print('Ephemeris was split in {} parts for better search of stars'.format(len(divisions)))

    # makes predictions for each division
    occs = []
    for i, vals in enumerate(divisions):
        nt = t[vals]
        if log:
            print('\nSearching occultations in part {}/{}'.format(i+1, len(divisions)))
            print("Generating Ephemeris between {} and {} ...".format(nt.min(), nt.max()))
        ncoord = ephem.get_position(nt)
        ra = np.mean([ncoord.ra.min().deg, ncoord.ra.max().deg])
        dec = np.mean([ncoord.dec.min().deg, ncoord.dec.max().deg])
        mindist = np.arcsin(radius/ncoord.distance).max() + sigma*np.max([ephem.error_ra.value, ephem.error_dec.value])*u.arcsec
        width = ncoord.ra.max() - ncoord.ra.min() + 2*mindist
        height = ncoord.dec.max() - ncoord.dec.min() + 2*mindist
        pos_search = SkyCoord(ra*u.deg, dec*u.deg)

        if log:
            print('Downloading stars ...')
        catalogue = vquery.query_region(pos_search, width=width, height=height, catalog='I/345/gaia2', cache=False)
        if len(catalogue) == 0:
            print('    No star found. The region is too small or VizieR is out.')
            continue
        catalogue = catalogue[0]
        if log:
            print('    {} Gaia-DR2 stars downloaded'.format(len(catalogue)))
            print('Identifying occultations ...')
        stars = SkyCoord(catalogue['RA_ICRS'], catalogue['DE_ICRS'])
        idx, d2d, d3d = stars.match_to_catalog_sky(ncoord)

        dist = np.arcsin(radius/ncoord[idx].distance) + sigma*np.max([ephem.error_ra.value, ephem.error_dec.value])*u.arcsec
        k = np.where(d2d < dist)[0]
        for ev in k:
            star = Star(code=catalogue['Source'][ev], nomad=False, log=False)
            c = star.geocentric(nt[idx][ev])
            pars = [star.code, SkyCoord(c.ra, c.dec), star.mag['G']]
            try:
                pars = np.hstack((pars, occ_params(star, ephem, nt[idx][ev])))
                occs.append(pars)
            except:
                pass

    meta = {'name': ephem.name, 'time_beg': time_beg, 'time_end': time_end, 'maglim': mag_lim, 'max_ca': mindist,
            'radius': ephem.radius.to(u.km).value, 'error_ra': ephem.error_ra.to(u.mas).value,
            'error_dec': ephem.error_dec.to(u.mas).value, 'ephem': ephem.meta['kernels']}
    if not occs:
        print('No stellar occultation was found.')
        return PredictionTable(meta=meta)
    # create astropy table with the params
    occs2 = np.transpose(occs)
    time = Time(occs2[3])
    geocentric = SkyCoord([ephem.get_position(time)])
    k = np.argsort(time)
    t = PredictionTable(
        time=time[k], coord_star=occs2[1][k], coord_obj=geocentric[k], ca=[i.value for i in occs2[4][k]],
        pa=[i.value for i in occs2[5][k]], vel=[i.value for i in occs2[6][k]], mag=occs2[2][k],
        dist=[i.value for i in occs2[7][k]], source=occs2[0][k], meta=meta)
    if log:
        print('{} occultations found.'.format(len(t)))
    return t


def xy2latlon(x, y, loncen, latcen, time):
    """ Calculates the longitude and latitude given projected positions x and y

    Parameters:
        x (int, float): Projected position in x, in the GCRS (meters)
        y (int, float): Projected position in y, in the GCRS (meters)
        loncen (int, float): Center longitude of projection (deg)
        latcen (int, float): Center latitude of projection (deg)
        time (Time): Time of refered projection

    Returns:
        lon, lat (float): Longitude and Latitude whose projection at loncen, lat results in x, y. (deg)
    """
    r = const.R_earth.to(u.m).value
    site_cen = EarthLocation(loncen*u.deg, latcen*u.deg)
    itrs_cen = site_cen.get_itrs(obstime=time)
    gcrs_cen = itrs_cen.transform_to(GCRS(obstime=time))

    z = np.array(y, ndmin=1)
    y = np.array(x, ndmin=1)
    x2 = r*r-y*y-z*z

    a = np.where(x2 >= 0.0)
    x = np.sqrt(x2[a])

    y = y[a]
    z = z[a]
    lon = np.repeat(1e+31, len(x2))
    lat = np.repeat(1e+31, len(x2))
    center_frame = SkyOffsetFrame(origin=gcrs_cen)
    if len(x) > 0:
        n = 0
        if not time.isscalar and len(time) == len(x2):
            time = time[a]
        while True:
            n += 1
            new_pos = SkyCoord(x*u.m, y*u.m, z*u.m, representation_type='cartesian', frame=center_frame[a])
            n_coord = new_pos.transform_to(GCRS(obstime=time))
            n_itrs = n_coord.transform_to(ITRS(obstime=time))
            n_site = n_itrs.earth_location
            n_site = EarthLocation(n_site.lon, n_site.lat, 0)
            itrs_site = n_site.get_itrs(obstime=time)
            gcrs_site = itrs_site.transform_to(GCRS(obstime=time))
            target1 = gcrs_site.transform_to(center_frame[a])
            if n == 4:
                lon[a] = n_site.lon.deg
                lat[a] = n_site.lat.deg
                break
            x = target1.cartesian.x.to(u.m).value
    return lon, lat


def latlon2xy(lon, lat, loncen, latcen):
    """ Calculates the projection of longitude and latitude in the loncen, latcen direction.

    Parameters:
        lon (int, float): Longitude to calculate projection
        lat (int, float): Latitude to calculate projection
        loncen (int, float): Center longitude of projection (deg)
        latcen (int, float): Center latitude of projection (deg)

    Returns:
        x, y (float): Projection of lon, lat at loncen, latcen, in the ITRS (meters)
    """
    site_cen = EarthLocation(loncen*u.deg, latcen*u.deg)
    itrs_cen = site_cen.get_itrs()

    lon = np.array(lon, ndmin=1)
    lat = np.array(lat, ndmin=1)
    site = EarthLocation(lon*u.deg, lat*u.deg, height=0*u.m)
    itrs_site = site.get_itrs()

    target = itrs_site.transform_to(SkyOffsetFrame(origin=itrs_cen))
    y = target.cartesian.y.to(u.m).value
    z = target.cartesian.z.to(u.m).value
    k = np.where(target.cartesian.x.to(u.m).value < 0.0)
    y[k] = 1e+31
    z[k] = 1e+31
    return y, z


def plot_occ_map(name, radius, coord, time, ca, pa, vel, dist, mag=0, longi=0, **kwargs):
    """ Plots map of the occultation

    Parameters:
        Required params:
        name (str): Name of the object
        radius (int, float): radius of the object (km)
        coord (str, SkyCoord): Coordinate of the star
            ("hh mm ss.sss dd mm ss.sss" or "hh.hhhhhhhh dd.dddddddd")
        time (str, Time): Instant of Closest Approach (iso or isot format)
        ca (int, float): Closest Approach Distance (arcsec)
        pa (int, float): Position Angle at C/A (deg)
        vel (int, vel): Velocity of the event (km/s)
        dist (int, float): Object distance at C/A (AU)

        Not required params (only printed in label):
        mag (int,float): Mag* = Normalized magnitude to vel=20km/s
        longi (int,float): East longitude of sub-planet point, deg, positive towards East

        Map configuration:
        nameimg (str): Change the name of the imaged saved.
        path (str): Path to a directory where to save map.
        resolution (int): Cartopy feature resolution. "1" means a resolution of "10m",
            "2" a resolution of "50m" and "3" a resolution of "100m". Default = 2
        states (bool): True to plot the state division of the countries. The states of
            some countries will only be shown depending on the resolution.
        zoom (int, float): Zooms in or out of the map.
        centermap_geo (list): Center the map given coordinates in longitude and latitude.
            It must be a list with two numbers. Default=None.
        centermap_delta (list): Displace the center of the map given displacement
            in X and Y, in km. It must be a list with two numbers. Default=None.
        centerproj (list): Rotates the Earth to show occultation with the center
            projected at a given longitude and latitude. It must be a list with two numbers
        labels (bool): Plots text above and below the map with the occultation parameters.
            Default=True.
        meridians (int): Plots lines representing the meridians for given interval. Default=30 deg
        parallels (int): Plots lines representing the parallels for given interval. Default=30 deg
        sites (dict): Plots site positions in map. It must be a python dictionary where the key is
            the name of the site, and the value is a list with longitude, latitude, delta_x,
            delta_y and color. delta_x and delta_y are displacement, in km, from the point
            of the site in the map and the name. color is the color of the point.
        site_name (bool): If True, it prints the name of the sites given, else it plots only the points
        countries (dict): Plots the names of countries. It must be a python dictionary where the key
            is the name of the country and the value is a list with longitude and latitude
            of the lower left part of the text.
        offset (list): applies an offset to the ephemeris, calculating new CA and instant of CA.
            It is a pair of delta_RA*cosDEC and delta_DEC.
        mapstyle (int): Define the color style of the map. 1 is the default black and white scale.
            2 is a colored map.
        error (int,float): Ephemeris error in mas. It plots a dashed line representing radius + error.
        ercolor (str): Changes the color of the lines of the error bar.
        ring (int,float): plots a dashed line representing the location of a ring.
            It is given in km, from the center.
        rncolor (str): Changes the color of ring lines.
        atm (int,float): plots a dashed line representing the location of an atmosphere.
            It is given in km, from the center.
        atcolor (str): Changes the color of atm lines.
        chord_delta (list): list with distances from center to plot chords
        chord_geo (2d-list): list with pairs of coordinates to plot chords
        chcolor (str): color of the line of the chords. Default: grey
        heights (list): It plots a circular dashed line showing the locations where the observer
            would observe the occultation at a given height above the horizons.
            This must be a list.
        hcolor (str): Changes the color of the height lines.
        mapsize (list): The size of figure, in cm. It must be a list with two values.
            Default = [46.0, 38.0].
        cpoints (int,float): Interval for the small points marking the center of shadow,
            in seconds. Default=60.
        ptcolor (str): Change the color of the center points.
        alpha (float): The transparency of the night shade, where 0.0 is full transparency
            and 1.0 is full black. Default = 0.2.
        fmt (str): The format to save the image. It is parsed directly by matplotlib.pyplot.
            Default = 'png'
        dpi (int): "Dots per inch". It defines the quality of the image. Default = 100.
        lncolor (str): Changes the color of the line that represents the limits of the shadow over Earth.
        outcolor (str): Changes the color of the lines that represents the limits of the shadow outside Earth
        nscale (int,float): Arbitrary scale for the size of the name of the site.
        cscale (int,float): Arbitrary scale for the name of the country.
        sscale (int,float): Arbitrary scale for the size of point of the site.
        pscale (int,float): Arbitrary scale for the size of the points that represent the center of the shadow
        arrow (bool): If true, it plots the arrow with the occultation direction.

        Comment: Only one of centermap_geo and centermap_delta can be given
    """
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature

    allowed_kwargs = ['alpha', 'arrow', 'atcolor', 'atm', 'centermap_delta', 'centermap_geo', 'centerproj',
                      'chcolor', 'chord_delta', 'chord_geo', 'countries', 'cpoints', 'cscale', 'dpi', 'ercolor',
                      'error', 'fmt', 'hcolor', 'heights', 'labels', 'lncolor', 'mapsize', 'mapstyle', 'meridians',
                      'nameimg', 'nscale', 'offset', 'outcolor', 'parallels', 'path', 'pscale', 'ptcolor',
                      'resolution', 'ring', 'rncolor', 'site_name', 'sites', 'sscale', 'states', 'zoom']
    input_tests.check_kwargs(kwargs, allowed_kwargs=allowed_kwargs)

    if not type(name) == str:
        raise TypeError('name keyword must be a string')

    radius = radius*u.km

    occs = {}
    try:
        occs['stars'] = SkyCoord(coord, frame='icrs', unit=(u.hourangle, u.degree))
    except:
        raise KeyError('"star" keyword is not in the format: "hh mm ss.sss dd mm ss.sss" or "hh.hhhhhhhh dd.dddddddd"')

    try:
        occs['datas'] = Time(time)
    except:
        raise KeyError('"time" keyword is not a iso or isot time format')
    occs['ca'] = ca*u.arcsec
    occs['posa'] = pa*u.deg
    occs['vel'] = vel*(u.km/u.s)
    occs['dist'] = dist*u.AU
    occs['magG'] = mag
    occs['longi'] = longi

    mapstyle = kwargs.get('mapstyle', 1)
    if mapstyle not in [1, 2]:
        raise ValueError('mapstyle must be 1 or 2]')

    resolution = kwargs.get('resolution', 2)
    if resolution not in [1, 2, 3]:
        raise TypeError('resolution keyword must be one of these: [1, 2, 3] where 1=10m, 2=50m and 3=100m')
    res = ['10m', '50m', '110m']
    resolution = res[resolution-1]

    nameimg = kwargs.get('nameimg', '{}_{}'.format(name, occs['datas'].isot))
    fmt = kwargs.get('fmt', 'png')
    dpi = kwargs.get('dpi', 100)
    step = kwargs.get('step', 1)
    mapsize = kwargs.get('mapsize', [46.0, 38.0])*u.cm
    erro = kwargs.get('error', None)
    ring = kwargs.get('ring', None)
    atm = kwargs.get('atm', None)
    cpoints = kwargs.get('cpoints', 60)
    states = kwargs.get('states', True)
    labels = kwargs.get('labels', True)
    meridians = kwargs.get('meridians', 30)
    parallels = kwargs.get('parallels', 30)
    nscale = kwargs.get('nscale', 1)
    cscale = kwargs.get('cscale', 1)
    sscale = kwargs.get('sscale', 1)
    pscale = kwargs.get('pscale', 1)
    heights = np.array(kwargs.get('heights'), None)
    alpha = kwargs.get('alpha', 0.2)
    centermap_geo = kwargs.get('centermap_geo', None)
    centermap_delta = kwargs.get('centermap_delta', None)
    if 'centermap_geo' in kwargs and 'centermap_delta' in kwargs:
        raise ValueError('User must give "centermap_geo" OR "centermap_delta"')
    zoom = kwargs.get('zoom', 1)
    off_ra, off_de = kwargs.get('offset', [0.0, 0.0])*u.mas
    arrow = kwargs.get('arrow', True)
    site_name = kwargs.get('site_name', True)
    path = kwargs.get('path', '.')
    chord_delta = np.array(kwargs.get('chord_delta', []), ndmin=1)*u.km
    chord_geo = kwargs.get('chord_geo', [])
    if len(chord_geo) > 0:
        try:
            b = np.array(chord_geo, ndmin=2)
            chord_geo = b.reshape(len(b), 2)
        except:
            raise ValueError('chord_geo must a set of pairs with longitude and latitude')
        chord_geo = EarthLocation(*chord_geo.T)

    sites = {}
    if 'sites' in kwargs.keys():
        if type(kwargs['sites']) == str and os.path.isfile(kwargs['sites']):
            data = np.loadtxt(kwargs['sites'], dtype={'names': ('name', 'lon', 'lat', 'offx', 'offy', 'color'),
                                                      'formats': ('S30', 'f8', 'f8', 'f8', 'f8', 'S30')},
                              delimiter=',', ndmin=1)
            for i, s in enumerate(data):
                sites[s['name'].strip().decode()] = [s['lon'], s['lat'], s['offx'], s['offy'], s['color'].strip().decode()]
        elif type(kwargs['sites']) == dict:
            sites = kwargs['sites']
        else:
            raise TypeError('sites keyword must be a file or a dictionary')

    countries = {}
    if 'countries' in kwargs.keys():
        if type(kwargs['countries']) == str and os.path.isfile(kwargs['countries']):
            data = np.loadtxt(kwargs['countries'], dtype={'names': ('name', 'lon', 'lat'), 'formats': ('S30', 'f8', 'f8')},
                              delimiter=',', ndmin=1)
            for i, c in enumerate(data):
                countries[c['name'].strip().decode()] = [c['lon'], c['lat']]
        elif type(kwargs['countries']) == dict:
            countries = kwargs['countries']
        else:
            raise TypeError('country keyword must be a file or a dictionary')

# calculates offsets
    dca = off_ra*np.sin(occs['posa']) + off_de*np.cos(occs['posa'])
    dt = (-(off_ra*np.cos(occs['posa']) - off_de*np.sin(occs['posa'])).to(u.rad)*occs['dist'].to(u.km)/np.absolute(occs['vel'])).value*u.s
    ca1 = occs['ca'] + dca
    data = occs['datas'] + dt

# define map parameters
    center_gcrs = GCRS(occs['stars'].ra, occs['stars'].dec, 1*u.R_earth, obstime=data)
    center_itrs = center_gcrs.transform_to(ITRS(obstime=data))
    center_map = center_itrs.earth_location
    centert = True
    if 'centerproj' in kwargs.keys():
        if type(kwargs['centerproj']) == EarthLocation:
            center_map = kwargs['centerproj']
        elif np.array(kwargs['centerproj']).shape == (2,):
            center_map = EarthLocation.from_geodetic(*kwargs['centerproj'], 0.0)
        else:
            raise TypeError('centerproj must be an Astropy EarthLocation Object or an array with Longitude and Latitude only')
        centert = False
    fig = plt.figure(figsize=(mapsize.to(u.imperial.inch).value))
    projection = ccrs.Orthographic(central_longitude=center_map.lon.value, central_latitude=center_map.lat.value)
    if labels:
        axf = plt.axes(projection=projection)
    else:
        axf = plt.axes([-0.001, -0.001, 1.002, 1.002], projection=projection)
    axf.set_global()

# calculates regions for zoom
    limits = None
    r = const.R_earth.to(u.m).value
    if centermap_geo is not None:
        cx, cy = latlon2xy(centermap_geo[0], centermap_geo[1], center_map.lon.value, center_map.lat.value)
        limits = [cx/1000.0, cy/1000.0]
        if np.any(np.absolute(limits) > r):
            raise ValueError('Coordinates for centermap_geo are outside the visible range.')
    elif centermap_delta is not None:
        limits = centermap_delta
    elif zoom != 1:
        limits = [0, 0]
    if limits is not None:
        dr = r/zoom
        l0 = (limits[0]*u.km).to(u.m).value
        l1 = (limits[1]*u.km).to(u.m).value
        dmsize = mapsize[0]/mapsize[1]
        if mapsize[1] < mapsize[0]:
            lx = l0 - dr*dmsize
            ux = l0 + dr*dmsize
            ly = l1 - dr
            uy = l1 + dr
        else:
            lx = l0 - dr
            ux = l0 + dr
            ly = l1 - dr/dmsize
            uy = l1 + dr/dmsize
        axf.set_xlim(lx, ux)
        axf.set_ylim(ly, uy)
        if labels and zoom > 1:
            centert = False

# plots features
    axf.coastlines(resolution=resolution, color='0.3')
    ocean = cfeature.NaturalEarthFeature('physical', 'ocean', resolution)
    land = cfeature.NaturalEarthFeature('physical', 'land', resolution)
    border = cfeature.NaturalEarthFeature('cultural', 'admin_0_countries', resolution)
    if mapstyle == 1:
        axf.add_feature(ocean, zorder=0, color='0.9')
        axf.add_feature(land, zorder=0, edgecolor='None', color='1.0')
        axf.add_feature(border, zorder=0.1, edgecolor='0.4', facecolor='None')
        axf.add_feature(cfeature.RIVERS, zorder=0, edgecolor='0.7')
        axf.add_feature(cfeature.LAKES, zorder=0, color='0.7')
        ptcolor = 'black'
        lncolor = 'blue'
        ercolor = 'blue'
        rncolor = 'blue'
        atcolor = 'blue'
        outcolor = 'red'
        hcolor = 'black'
        chcolor = 'gray'
    elif mapstyle == 2:
        axf.add_feature(ocean, zorder=0, facecolor=cfeature.COLORS['water'])
        axf.add_feature(land, zorder=0, edgecolor='None', facecolor=cfeature.COLORS['land'])
        axf.add_feature(border, zorder=0, edgecolor='0.5', facecolor=cfeature.COLORS['land'])
        axf.add_feature(border, zorder=0.1, edgecolor='0.5', facecolor='None')
        axf.add_feature(cfeature.RIVERS, zorder=0)
        axf.add_feature(cfeature.LAKES, zorder=0)
        ptcolor = 'red'
        lncolor = 'blue'
        ercolor = 'red'
        rncolor = 'black'
        atcolor = 'black'
        outcolor = 'red'
        hcolor = 'black'
        chcolor = 'gray'
    if states:
        states_r = cfeature.NaturalEarthFeature('cultural', 'admin_1_states_provinces', resolution)
        axf.add_feature(states_r, zorder=0, edgecolor='0.6', facecolor='None')

    gl = axf.gridlines(xlocs=np.arange(-180, 180.001, meridians), ylocs=np.arange(-90, 90.001, parallels))
    gl.n_steps = 180
    sun = get_sun(data)
    sun_lat = sun.dec
    sun_lon = sun.ra - data.sidereal_time('mean', 'greenwich')
    pole_lon = sun_lon.deg
    pole_lat = sun_lat.deg
    proj_sun = ccrs.Orthographic(central_longitude=pole_lon+180, central_latitude=-pole_lat)
    bordx = r*np.cos(np.arange(0, 361, 0.5)*u.deg)
    bordy = r*np.sin(np.arange(0, 361, 0.5)*u.deg)
    axf.fill(bordx, bordy, transform=proj_sun, linewidth=0, color='black', alpha=alpha)
    axf.fill(bordx*np.cos(18*u.deg), bordy*np.cos(18*u.deg), transform=proj_sun, linewidth=0, color='black', alpha=alpha)

    ptcolor = kwargs.get('ptcolor', ptcolor)
    lncolor = kwargs.get('lncolor', lncolor)
    ercolor = kwargs.get('ercolor', ercolor)
    rncolor = kwargs.get('rncolor', rncolor)
    atcolor = kwargs.get('atcolor', atcolor)
    outcolor = kwargs.get('outcolor', outcolor)
    hcolor = kwargs.get('hcolor', hcolor)
    chcolor = kwargs.get('chcolor', chcolor)

# calculates path
    vec = np.arange(0, int(8000/(np.absolute(occs['vel'].value))), step)
    vec = np.sort(np.concatenate((vec, -vec[1:]), axis=0))
    pa = Angle(occs['posa'])
    pa.wrap_at('180d', inplace=True)
    if pa > 90*u.deg:
        paplus = pa - 180*u.deg
    elif pa < -90*u.deg:
        paplus = pa + 180*u.deg
    else:
        paplus = pa
    deltatime = vec*u.s
    datas1 = data + deltatime
    centers_gcrs = GCRS(np.repeat(occs['stars'].ra, len(datas1)), np.repeat(occs['stars'].dec, len(datas1)),
                        1*u.R_earth, obstime=datas1)
    centers_itrs = centers_gcrs.transform_to(ITRS(obstime=datas1))
    centers = centers_itrs.earth_location

    dista = (occs['dist'].to(u.km)*ca1.to(u.rad)).value*u.km
    ax = dista*np.sin(pa) + (deltatime*occs['vel'])*np.cos(paplus)
    by = dista*np.cos(pa) - (deltatime*occs['vel'])*np.sin(paplus)

    ax2 = ax - (radius)*np.sin(paplus)
    by2 = by - (radius)*np.cos(paplus)
    lon1, lat1 = xy2latlon(ax2.to(u.m).value, by2.to(u.m).value, centers.lon.value, centers.lat.value, datas1)
    j = np.where(lon1 < 1e+30)
    axf.plot(lon1[j], lat1[j], transform=ccrs.Geodetic(), color=lncolor)
    j = np.where(lon1 > 1e+30)
    if 'centerproj' not in kwargs:
        plt.plot(ax2[j].to(u.m).value, by2[j].to(u.m).value, color=outcolor, clip_on=(not centert), zorder=-0.2)

    ax3 = ax + (radius)*np.sin(paplus)
    by3 = by + (radius)*np.cos(paplus)
    lon2, lat2 = xy2latlon(ax3.to(u.m).value, by3.to(u.m).value, centers.lon.value, centers.lat.value, datas1)
    j = np.where(lon2 < 1e+30)
    axf.plot(lon2[j], lat2[j], transform=ccrs.Geodetic(),  color=lncolor)
    j = np.where(lon2 > 1e+30)
    if 'centerproj' not in kwargs:
        plt.plot(ax3[j].to(u.m).value, by3[j].to(u.m).value, color=outcolor, clip_on=(not centert), zorder=-0.2)

# plots chords_delta
    for val in chord_delta:
        ax2 = ax + val*np.sin(paplus)
        by2 = by + val*np.cos(paplus)
        lon1, lat1 = xy2latlon(ax2.to(u.m).value, by2.to(u.m).value, centers.lon.value, centers.lat.value, datas1)
        j = np.where(lon1 < 1e+30)
        axf.plot(lon1[j], lat1[j], transform=ccrs.Geodetic(), color=chcolor)

# plots chords_geo
    for coord_geo in chord_geo:
        xt, yt = latlon2xy(coord_geo.lon.deg, coord_geo.lat.deg, centers.lon.value, centers.lat.value)*u.m
        val = np.sqrt((xt-ax)**2 + (yt-by)**2)
        k = val.argmin()
        ang = np.arctan2((yt-by)[k], (xt-ax)[k])
        val = np.sign(np.sin(ang))*val[k]
        ax2 = ax + val*np.sin(paplus)
        by2 = by + val*np.cos(paplus)
        lon1, lat1 = xy2latlon(ax2.to(u.m).value, by2.to(u.m).value, centers.lon.value, centers.lat.value, datas1)
        j = np.where(lon1 < 1e+30)
        axf.plot(lon1[j], lat1[j], transform=ccrs.Geodetic(), color=chcolor)

# plots error
    if erro is not None:
        err = erro*u.mas
        errd = (occs['dist'].to(u.km)*err.to(u.rad)).value*u.km
        ax2 = ax - errd*np.sin(paplus) - radius*np.sin(paplus)
        by2 = by - errd*np.cos(paplus) - radius*np.cos(paplus)
        lon1, lat1 = xy2latlon(ax2.to(u.m).value, by2.to(u.m).value, centers.lon.value, centers.lat.value, datas1)
        j = np.where(lon1 < 1e+30)
        axf.plot(lon1[j], lat1[j], '--', transform=ccrs.Geodetic(),  color=ercolor)

        ax3 = ax + errd*np.sin(paplus) + radius*np.sin(paplus)
        by3 = by + errd*np.cos(paplus) + radius*np.cos(paplus)
        lon2, lat2 = xy2latlon(ax3.to(u.m).value, by3.to(u.m).value, centers.lon.value, centers.lat.value, datas1)
        j = np.where(lon2 < 1e+30)
        axf.plot(lon2[j], lat2[j], '--', transform=ccrs.Geodetic(),  color=ercolor)

# plots ring
    if ring is not None:
        rng = ring*u.km
        ax2 = ax - rng*np.sin(paplus)
        by2 = by - rng*np.cos(paplus)
        lon1, lat1 = xy2latlon(ax2.to(u.m).value, by2.to(u.m).value, centers.lon.value, centers.lat.value, datas1)
        j = np.where(lon1 < 1e+30)
        axf.plot(lon1[j], lat1[j], '--', transform=ccrs.Geodetic(),  color=rncolor)

        ax3 = ax + rng*np.sin(paplus)
        by3 = by + rng*np.cos(paplus)
        lon2, lat2 = xy2latlon(ax3.to(u.m).value, by3.to(u.m).value, centers.lon.value, centers.lat.value, datas1)
        j = np.where(lon2 < 1e+30)
        axf.plot(lon2[j], lat2[j], '--', transform=ccrs.Geodetic(),  color=rncolor)

# plots atm
    if atm is not None:
        atmo = atm*u.km
        ax2 = ax - atmo*np.sin(paplus)
        by2 = by - atmo*np.cos(paplus)
        lon1, lat1 = xy2latlon(ax2.to(u.m).value, by2.to(u.m).value, centers.lon.value, centers.lat.value, datas1)
        j = np.where(lon1 < 1e+30)
        axf.plot(lon1[j], lat1[j], '--', transform=ccrs.Geodetic(),  color=atcolor)

        ax3 = ax + atmo*np.sin(paplus)
        by3 = by + atmo*np.cos(paplus)
        lon2, lat2 = xy2latlon(ax3.to(u.m).value, by3.to(u.m).value, centers.lon.value, centers.lat.value, datas1)
        j = np.where(lon2 < 1e+30)
        axf.plot(lon2[j], lat2[j], '--', transform=ccrs.Geodetic(),  color=atcolor)

# plots center points
    vec = np.arange(0, int(8000/(np.absolute(occs['vel'].value))), cpoints)
    deltatime = np.sort(np.concatenate((vec, -vec[1:]), axis=0))*u.s
    axc = dista*np.sin(pa) + (deltatime*occs['vel'])*np.cos(paplus)
    byc = dista*np.cos(pa) - (deltatime*occs['vel'])*np.sin(paplus)
    plt.plot(axc.to(u.m).value, byc.to(u.m).value, 'o', color=ptcolor, clip_on=(not centert),
             markersize=mapsize[0].value*pscale*8.0/46.0, zorder=-0.2)

    datas2 = data + deltatime
    centers_p_gcrs = GCRS(np.repeat(occs['stars'].ra, len(datas2)), np.repeat(occs['stars'].dec, len(datas2)),
                          1*u.R_earth, obstime=datas2)
    centers_p_itrs = centers_p_gcrs.transform_to(ITRS(obstime=datas2))
    centers_p = centers_p_itrs.earth_location
    clon1, clat1 = xy2latlon(axc.to(u.m).value, byc.to(u.m).value, centers_p.lon.value, centers_p.lat.value, datas2)
    j = np.where(clon1 < 1e+30)
    axf.plot(clon1[j], clat1[j], 'o', transform=ccrs.Geodetic(), color=ptcolor, clip_on=True,
             markersize=mapsize[0].value*pscale*8.0/46.0)

    datas1 = data + deltatime
    center_gcrs = GCRS(np.repeat(occs['stars'].ra, 1), np.repeat(occs['stars'].dec, 1),
                       1*u.R_earth, obstime=data)
    center_itrs = center_gcrs.transform_to(ITRS(obstime=data))
    center = center_itrs.earth_location
    xp = [(dista.to(u.m)*np.sin(pa)).value]
    yp = [(dista.to(u.m)*np.cos(pa)).value]
    loncen, latcen = xy2latlon(xp, yp, center.lon.value, center.lat.value, data)
    j = np.where(loncen < 1e+30)
    if len(j) > 0:
        axf.plot(loncen[j], latcen[j], 'o', transform=ccrs.Geodetic(), color=ptcolor, clip_on=True,
                 markersize=mapsize[0].value*pscale*24.0/46.0)
    elif not centert:
        plt.plot(xp, yp, 'o', color=ptcolor, clip_on=False, markersize=mapsize[0].value*pscale*24.0/46.0)

# plots the heights
    if 'heights' in kwargs.keys():
        for h in heights:
            plt.plot(bordx*np.cos(h*u.deg), bordy*np.cos(h*u.deg), linestyle='dotted', color=hcolor)

# plots the the direction arrow
    if arrow:
        if limits is None:
            plt.quiver(5500000, -5500000, (np.sin(paplus+90*u.deg)*np.sign(occs['vel'])).value,
                       (np.cos(paplus+90*u.deg)*np.sign(occs['vel'])).value, width=0.005)
        else:
            plt.quiver(lx + (ux-lx)*0.9, ly + (uy-ly)*0.1, (np.sin(paplus+90*u.deg)*np.sign(occs['vel'])).value,
                       (np.cos(paplus+90*u.deg)*np.sign(occs['vel'])).value, width=0.005, zorder=1.3)

# plots the countries names
    for country in countries.keys():
        plt.text(countries[country][0], countries[country][1], country, transform=ccrs.Geodetic(),
                 weight='bold', color='grey', fontsize=30*cscale, family='monospace')

# plots the sites
    for site in sites.keys():
        s = EarthLocation.from_geodetic(sites[site][0], sites[site][1], 0.0*u.km)
        axf.plot(s.lon.deg, s.lat.deg, 'o', transform=ccrs.Geodetic(),
                 markersize=mapsize[0].value*sscale*10.0/46.0, color=sites[site][4])
        if site_name:
            xt, yt = latlon2xy(s.lon.deg, s.lat.deg, center_map.lon.value, center_map.lat.value)
            axf.text(xt + sites[site][2]*1000, yt+sites[site][3]*1000, site, weight='bold',
                     fontsize=25*nscale, family='monospace')

# Define the title and label of the output
    title = ('Object        Diam   Tmax   dots <> ra_offset_dec\n'
             '{:10s} {:4.0f} km  {:5.1f}s  {:02d} s <>{:+6.1f} {:+6.1f} \n'.
             format(name, 2*radius.value, (2*radius/np.absolute(occs['vel'])).value,
                    cpoints, off_ra.value, off_de.value))
    labelx = ("\n year-m-d    h:m:s UT     ra__dec__J2000__candidate    C/A    P/A    vel   Delta   G*  long\n"
              "{}  {:02d} {:02d} {:07.4f} {:+03d} {:02d} {:06.3f} {:6.3f} {:6.2f} {:6.2f}  {:5.2f} {:5.1f}  {:3.0f}".
              format(data.iso, int(occs['stars'].ra.hms.h), int(occs['stars'].ra.hms.m), occs['stars'].ra.hms.s,
                     int(occs['stars'].dec.dms.d), np.absolute(int(occs['stars'].dec.dms.m)),
                     np.absolute(occs['stars'].dec.dms.s), ca1.value, occs['posa'].value,
                     occs['vel'].value, occs['dist'].value, occs['magG'], occs['longi']))

# plots the map
    if labels:
        axf.set_title(title, family='monospace', weight='bold', fontsize=22)
        axf.text(0.5, -0.1, labelx, va='bottom', ha='center', rotation='horizontal', rotation_mode='anchor',
                 transform=axf.transAxes, family='monospace', weight='bold', fontsize=22)
    filepath = os.path.join(path, '{}.{}'.format(nameimg, fmt))
    plt.savefig(filepath, format=fmt, dpi=dpi)
    print('{}.{} generated'.format(nameimg, fmt))
    plt.clf()
    plt.close()
