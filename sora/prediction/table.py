import glob
import os
import warnings

import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord, get_sun, get_moon
from astropy.table import Table, Row, Column
from astropy.time import Time

from sora.config import input_tests

__all__ = ['PredictionTable']


class PredictRow(Row):
    """An Astropy Row object modified for Prediction purposes.
    """

    def plot_occ_map(self, **kwargs):
        """
        Parameters
        ----------
        radius : `int`, `float`
            The radius of the shadow. If not given it uses saved value.

        nameimg : `str`
            Change the name of the image saved.

        path : `str`
            Path to a directory where to save map.

        resolution : `int`, default=2
        Cartopy feature resolution.\n
        - ``1`` means a resolution of "10m";\n
        - ``2`` a resolution of "50m";\n
        - ``3`` a resolution of "100m".

        states : `bool`
            If True, plots the states borders of the countries. The states
            of some countries will only be shown depending on the resolution.

        zoom : `int`, `float`
            Zooms in or out of the map.

        centermap_geo : `list`, default=None
            Center the map given coordinates in longitude and latitude. It must be
            a list with two numbers.

        centermap_delta : `list`, default=None
            Displace the center of the map given displacement in X and Y, in km.
            It must be a list with two numbers.

        centerproj : `list`
            Rotates the Earth to show occultation with the center projected at a
            given longitude and latitude. It must be a list with two numbers.

        labels : `bool`, default=True
            Plots text above and below the map with the occultation parameters.

        meridians : `int`, default=30
            Plots lines representing the meridians for given interval, in degrees.

        parallels : `int`, default=30
            Plots lines representing the parallels for given interval, in degrees.

        sites : `dict`
            Plots site positions in map. It must be a python dictionary where the
            key is  the `name` of the site, and the value is a list with `longitude`,
            `latitude`, `delta_x`, `delta_y` and `color`. `delta_x` and `delta_y`
            are displacement, in km, from the point position of the site in the map
            and the `name`. `color` is the color of the point.

        site_name : `bool`
            If True, it prints the name of the sites given, else it plots only the points.

        countries : `dict`
            Plots the names of countries. It must be a python dictionary where the
            key is the name of the country and the value is a list with longitude
            and latitude of the lower left part of the text.

        offset : `list`
            Applies an offset to the ephemeris, calculating new CA and instant of
            CA. It is a pair of `delta_RA*cosDEC` and `delta_DEC`.

        mapstyle : `int`, default=1
            Define the color style of the map. ``'1'`` is the default black
            and white scale. ``'2'`` is a colored map.

        error : `int`, `float`
            Ephemeris error in mas. It plots a dashed line representing radius + error.

        ercolor : `str`
            Changes the color of the lines of the error bar.

        ring : `int`, `float`
            Plots a dashed line representing the location of a ring. It is given
            in km, from the center.

        rncolor : `str`
            Changes the color of ring lines.

        atm : `int`, `float`
            Plots a dashed line representing the location of an atmosphere. It is
            given in km, from the center.

        atcolor : `str`
            Changes the color of atm lines.

        chord_delta : `list`
            List with distances from center to plot chords.

        chord_geo : `2d-list`
            List with pairs of coordinates to plot chords.

        chcolor : `str`, default='grey'
            Color of the line of the chords.

        heights : `list`
            It plots a circular dashed line showing the locations where the observer
            would observe the occultation at a given height above the horizons.
            This must be a list.

        hcolor : `str`
            Changes the color of the height lines.

        mapsize : `list`, default= [46.0, 38.0]
            The size of figure, in cm. It must be a list with two values.

        cpoints : `int`, `float`, default=60
            Interval for the small points marking the center of shadow, in seconds.

        ptcolor : `str`
            Change the color of the center points.

        alpha : `float`, default=0.2
            The transparency of the night shade, where 0.0 is full transparency and
            1.0 is full black.

        fmt : `str`, default:'png'
            The format to save the image. It is parsed directly by `matplotlib.pyplot`.

        dpi : `int`, default=100
            Resolution in "dots per inch". It defines the quality of the image.

        lncolor : `str`
            Changes the color of the line that represents the limits of the shadow
            over Earth.

        outcolor :`str`
            Changes the color of the lines that represents the limits of the shadow
            outside Earth.

        scale : `int`, `float`
            Arbitrary scale for the size of the name of the site.

        cscale : `int`, `float`
            Arbitrary scale for the name of the country.

        sscale : `int`, `float`
            Arbitrary scale for the size of point of the site.

        pscale : `int`, `float`
            Arbitrary scale for the size of the points that represent the center of
            the shadow.

        arrow : `bool`
            If True, it plots the arrow with the occultation direction.


        Note
        ----
        Only one of centermap_geo and centermap_delta can be given.
        """
        from .occmap import plot_occ_map

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
                values['G'] = Column(kwargs['mag'] * u.mag, format='6.3f', unit='mag')
                del kwargs['mag']
            else:
                values['G'] = kwargs['mag_20'] - 2.5*np.log10(np.absolute(values['Vel'])/20.0)
                values['G'].unit = 'mag'
                values['G'].format = '6.3f'
            if 'mag_20' in kwargs.keys():
                values['G*'] = Column(kwargs['mag_20'] * u.mag, format='6.3f', unit='mag')
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
            catalogue = kwargs.get('meta', {}).get('catalogue', '')
            if 'source' in kwargs.keys():
                values[f'{catalogue} Source ID'] = Column(kwargs['source'], dtype=np.int64)
                del kwargs['source']
            else:
                values[f'{catalogue} Source ID'] = Column(np.repeat('', len(time)))
            super().__init__(values, **kwargs)
        else:
            super().__init__(*args, **kwargs)

    def __itens_by_epoch(self, date):
        """ Gets item list for all occultations that matches the given date

        Parameters
        ----------
        date : `str`
            Date to match

        Returns
        -------
        item : `list`
            The list of occultations that matches the date.
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
        """Creates a PredictionTable Table reading from a PRAIA table.

        Parameters
        ----------
        filename : `str`
            Path to the PRAIA table file.

        name : `str`
            Name of the Object of the prediction.

        radius : `int`, `float`, optional
            Object radius, in km.
            If not given it's searched in online database.
            When not found online, the default is set to zero.

        Returns
        -------
         : `sora.prediction.PredictionTable`
            A PredictionTable object.
        """
        from sora.body.utils import search_satdb, search_sbdb
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
        try:
            data = search_satdb(name.lower())
            radius = data.get('diameter', 0) / 2
        except ValueError:
            try:
                data = search_sbdb(name.lower())
                radius = data.get('diameter', 0)/2
            except ValueError:
                radius = 0
        error_ra, error_dec = 0, 0
        radius = kwargs.get('radius', radius)*u.km
        meta = {'name': name, 'radius': radius, 'max_ca': max_ca, 'ephem': lines[17].split()[-1],
                'error_ra': error_ra*1000, 'error_dec': error_dec*1000}
        return cls(time=time, coord_star=coord_star, coord_obj=coord_obj, ca=dados['ca'], pa=dados['pa'],
                   vel=dados['vel'], mag_20=dados['mR'], dist=dados['delta'],
                   long=dados['long'], loct=dados['loct'], meta=meta)

    def to_praia(self, filename):
        """Writes PredictionTable to PRAIA format.

        Parameters
        ----------
        filename : `str`
            Name of the file to save table.
        """
        from .values import praia_occ_head
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
        """Writes PredictionTable to OccultWatcher feeder update file format.
        Tables will be saved in two files: "tableOccult_update.txt" and "LOG.dat"

        Parameters
        ----------
        ow_des : `str`
            Occult Watcher designation for the object.

        mode : `str`, default='append'
            Use ``'append'`` to append table to already existing file and
            ``'restart'`` to overwrite existing file.
        """
        from .values import ow_occ_head
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
        basename = kwargs.get('nameimg', None)
        for i in range(len(self)):
            if basename and len(self) > 1:
                kwargs['nameimg'] = f'{basename}_{i}'
            self[i].plot_occ_map(**kwargs)

    plot_occ_map.__doc__ = PredictRow.plot_occ_map.__doc__

    def remove_occ(self, date):
        """Removes stellar occultations from table.

        Parameters
        ----------
        date : `str`, `list`
            Date or list of dates of the occultation to be removed.
            The dates mut be as shown in the 'Epoch' column. If the date is not
            complete, the function will select all occultations that matches the
            given string. For instance, ``date='2020-06'`` will remove all
            occultations from the month of June 2020.
        """
        if type(date) == str:
            date = [date]
        itens = list(set([y for d in date for y in self.__itens_by_epoch(d)]))
        self.remove_rows(itens)

    def keep_from_selected_images(self, path='.'):
        """Keeps predictions which images were not deleted in given path.
        This function uses the name of the images to identify predictions.
        The name must be the automatic one generated by plot_occ_map().
        The format of the image is not relevant.

        Parameters
        ----------
        path : `str`
            Path where images are located.
        """
        itens = []
        for i, tca in enumerate(self['Epoch']):
            name = '{}_{}'.format(self.meta['name'], tca.isot)
            if not glob.glob(os.path.join(path, name)+'*'):
                itens.append(i)
        self.remove_rows(itens)
