import warnings

import astropy.units as u
import numpy as np
from astropy.coordinates import SkyOffsetFrame
from astropy.time import Time

__all__ = ['Chord']


class Chord:
    """Defines an Occultation Chord.

    Attributes
    ----------
    name : `str`
        The name of the Chord.

    observer : `sora.Observer`
        The site of observation.

    lightcurve : `sora.LightCurve`
        The lightcurve observed.

    """

    def __init__(self, *, name, observer, lightcurve):

        from sora.lightcurve import LightCurve
        from sora.observer import Observer, Spacecraft

        if not isinstance(observer, (Observer, Spacecraft)):
            raise ValueError('obs must be an Observer object')
        if not isinstance(lightcurve, LightCurve):
            raise ValueError('lightcurve must be a LightCurve object')
        self._shared_with = {}
        self._name = name
        self._observer = observer
        self._lightcurve = lightcurve
        self._isable = {}
        self._method = 'geocenter'

    @property
    def name(self):
        return self._name

    @property
    def observer(self):
        return self._observer

    @property
    def lightcurve(self):
        return self._lightcurve

    def status(self):
        """Returns if the chord is positive or negative
        """
        im = getattr(self.lightcurve, 'immersion', None)
        em = getattr(self.lightcurve, 'emersion', None)
        if im is None and em is None:
            return 'negative'
        else:
            return 'positive'

    @property
    def is_able(self):
        """Return a dictionary with the references enabled and disabled for the fit.
        """
        im = getattr(self.lightcurve, 'immersion', None)
        if im is None and 'immersion' in self._isable:
            del(self._isable['immersion'])
        if im is not None and 'immersion' not in self._isable:
            self._isable['immersion'] = True
        em = getattr(self.lightcurve, 'emersion', None)
        if em is None and 'emersion' in self._isable:
            del(self._isable['emersion'])
        if em is not None and 'emersion' not in self._isable:
            self._isable['emersion'] = True
        return self._isable

    def enable(self, *, time=None):
        """Enables a contact point of the curve to be used in the fit.

        Parameters
        ----------
        time : `None`, `str`
            If ``None``, it will enable all contact points.
            If ``'immersion'`` or ``'emersion'``, it will enable respective
            contact point.
        """
        if time not in [None, 'immersion', 'emersion']:
            raise ValueError('time can only be "immersion", "emersion" or None for both values')
        isable = self.is_able
        for key in isable:
            if time is None or time == key:
                isable[key] = True

    def disable(self, *, time=None):
        """Disables a contact point of the curve to be used in the fit.

        Parameters
        ----------
        time : `None`, `str`
            If None, it will disable all contact points.
            If ``'immersion'`` or ``'emersion'``, it will disable respective
            contact point.
        """
        if time not in [None, 'immersion', 'emersion']:
            raise ValueError('time can only be "immersion", "emersion" or None for both values')
        isable = self.is_able
        for key in isable:
            if time is None or time == key:
                isable[key] = False

    def get_fg(self, *, time=None, vel=False):
        """Returns the on-sky position for this chord at a given time.

        This Chord object must be associated to an Occultation to work, since
        it needs the position of the star and an ephemeris.

        Parameters
        ----------
        time : `astropy.time.Time`, `str`
            It must a time or a list of time to calculate the on-sky position.
            If time is ``'immersion'``, ``'emersion'``, ``'start'`` or ``'end'``,
            it will get the respective on-sky position. If ``time=None``, it
            will return the on-sky position for the immersion and emersion times
            if positive or raise an error if negative.

        vel : `bool`
            If True, it will return the on-sky velocity as well.

        Return
        ------
        f, g, [vf, vg] : `float`, `float`, [`float`, `float`]
            The geocentric on-sky Orthographic projection of the object with `f`
            pointing to the celestial North and `g` pointing to the celestial
            East. The respective velocities (vf, vg) are returned if ``vel=True``.
        """
        from sora.observer import Spacecraft

        if 'chordlist' not in self._shared_with:
            raise ValueError('{} must be associated to an Occultation to use this function'.format(self.__class__.__name__))
        occtime = self._shared_with['chordlist']['time']
        star = self._shared_with['chordlist']['star']
        coord = star.get_position(time=occtime, observer='geocenter')
        ephem = self._shared_with['chordlist']['ephem']

        ref_times = {'immersion': 'immersion', 'emersion': 'emersion', 'start': 'initial_time', 'end': 'end_time'}

        if isinstance(time, str) and time in ref_times.keys():
            time = getattr(self.lightcurve, ref_times[time])
        elif time is None:
            if self.status() == 'negative':
                raise ValueError("{} {}'s LightCurve does not have immersion or emersion time".format(
                    self.__class__.__name__, self.name))
            else:
                immersion = self.lightcurve.immersion
                emersion = self.lightcurve.emersion
                time = Time([immersion, emersion])
        else:
            time = Time(time)

        tca_diff = np.array(np.absolute((time - occtime).jd), ndmin=1)
        if any(tca_diff > 0.5):
            warnings.warn('The difference between a given time and the closest approach instant is {:.2f} days. '
                          'This position could not have a physical meaning.'.format(tca_diff.max()))

        if self._method == 'geocenter' and not isinstance(self.observer, Spacecraft):
            ksio1, etao1 = self.observer.get_ksi_eta(time=time, star=coord)
            ksie1, etae1 = ephem.get_ksi_eta(time=time, star=coord)
            f = ksio1 - ksie1
            g = etao1 - etae1
        else:
            pos_star = star.get_position(time=time, observer=self.observer)
            pos_ephem = ephem.get_position(time=time, observer=self.observer)
            _, f, g, = - pos_ephem.transform_to(SkyOffsetFrame(origin=pos_star)).cartesian.xyz.to(u.km).value

        if not vel:
            return f, g

        if self._method == 'geocenter' and not isinstance(self.observer, Spacecraft):
            ksio2, etao2 = self.observer.get_ksi_eta(time=time + 0.1 * u.s, star=coord)
            ksie2, etae2 = ephem.get_ksi_eta(time=time + 0.1 * u.s, star=coord)
            nf = ksio2 - ksie2
            ng = etao2 - etae2
        else:
            pos_star = star.get_position(time=time + 0.1 * u.s, observer=self.observer)
            pos_ephem = ephem.get_position(time=time + 0.1 * u.s, observer=self.observer)
            _, nf, ng, = - pos_ephem.transform_to(SkyOffsetFrame(origin=pos_star)).cartesian.xyz.to(u.km).value

        return f, g, (nf-f)/0.1, (ng-g)/0.1

    def path(self, *, segment='standard', step=1):
        """Returns the on-sky path of this chord.

        This Chord object must be associated to an Occultation to work, since it
        needs the position of the star and an ephemeris.

        Parameters
        ----------
        segment : `str`
            The segment to get the path of the chord. The available options are:

            ``'positive'``: to get the path between the immersion and emersion
            times if the chord is positive.

            ``'negative'``: to get the path between the start and end of
            observation if the chord is negative.

            ``'standard'``: to get the 'positive' path if the chord is positive
            or 'negative' if the chord is negative.

            ``'full'``: to get  the path between the start and end of observation
            independent if the chord is positive or negative.

            ``'outer'``: to get the path outside the 'positive' path, for
            instance between the start and immersion times and between the
            emersion and end times.

            ``'error'``: to get the path corresponding to the error bars. Be aware
            that some of these segments may return more than 1 path, for example
            ``segment='error'``.

        step : `int`, `float`, `str`
            If a number, it corresponds to the step, in seconds, for each point
            of the path.

            This correspond to an approximate value if it is not a multiple of
            the interval of the segment.

            The step can also be equal to ``'exposure'``. In this case, the path
            will return a set of pairs where each pair will have the position
            at the beginning of the exposure and the end of the exposure.


        Return
        ------
        path f, path g : `array`, `array`
            It will return the path separated by `f` and `g`. If ``step='exposure'``,
            it will return f and g intercalated for each exposure.
        """
        if 'chordlist' not in self._shared_with:
            raise ValueError('{} must be associated to an Occultation to use this function'.format(self.__class__.__name__))

        intervals = []
        try:
            immersion = self.lightcurve.immersion
            emersion = self.lightcurve.emersion
            if segment in ['positive', 'standard']:
                intervals.append((immersion, emersion))
            elif segment == 'outer':
                intervals.append((self.lightcurve.initial_time, immersion))
                intervals.append((emersion, self.lightcurve.end_time))
            elif segment == 'error':
                immersion_err = self.lightcurve.immersion_err*u.s
                emersion_err = self.lightcurve.emersion_err*u.s
                intervals.append((immersion-immersion_err, immersion+immersion_err))
                intervals.append((emersion-emersion_err, emersion+emersion_err))
            elif segment == 'full':
                intervals.append((self.lightcurve.initial_time, self.lightcurve.end_time))
            else:
                raise ValueError('segment {} can not be determined for chord {}'.format(segment, self.name))
        except AttributeError:
            if segment in ['negative', 'standard', 'full']:
                intervals.append((self.lightcurve.initial_time, self.lightcurve.end_time))
            else:
                raise ValueError('segment {} can not be determined for chord {}'.format(segment, self.name))
        exposure = False
        if step == 'exposure':
            exposure = True
            if segment == 'error':
                warnings.warn('"exposure" not available for segment "error".')
                exposure = False
                step = 1
            else:
                try:
                    exptime = self.lightcurve.exptime
                except AttributeError:
                    warnings.warn('{} {} does not have "exptime". Setting continuous path'.format(
                        self.__class__.__name__, self.name))
                    exposure = False
                    step = 1
        elif not isinstance(step, (int, float)):
            raise ValueError("step must be a number or the string 'exposure'")

        vals = []
        for interval in intervals:
            if exposure:
                tref = self.lightcurve.tref
                times = tref + self.lightcurve.time*u.s
                times = times[np.where((times >= interval[0]) & (times <= interval[1]))]
                time_beg = times - exptime/2.0*u.s
                time_end = times + exptime/2.0*u.s
                f1, g1 = self.get_fg(time=time_beg)
                f2, g2 = self.get_fg(time=time_end)
                v = np.array([f1, f2, g1, g2])
                nv = v.T.reshape(2*len(f1), 2)
                for i in nv:
                    vals.append(i)
            else:
                dt = interval[1] - interval[0]
                n = int(dt.sec/step)+1
                if n < 2:
                    n = 2
                time = interval[0] + np.linspace(0, dt.sec, n)*u.s
                f, g = self.get_fg(time=time)
                vals.append(f)
                vals.append(g)

        return vals

    def plot_chord(self, *, segment='standard', only_able=False, ax=None, linestyle='-', **kwargs):
        """Plots the on-sky path of this chord.

        This Chord object must be associated to an Occultation to work, since
        it needs the position of the star and an ephemeris.

        Parameters
        ----------
        segment : `str`
            The segment to plot the chord. The available options are:

            ``'positive'``: to get the path between the immersion and emersion
            times if the chord is positive.

            ``'negative'``: to get the path between the start and end of
            observation if the chord is negative.

            ``'standard'``: to get the 'positive' path if the chord is positive
                or 'negative' if the chord is negative.

            ``'full'``: to get  the path between the start and end of observation
            independent if the chord is positive or negative.

            ``'outer'``: to get the path outside the 'positive' path, for
            instance between the start and immersion times and between the
            emersion and end times.

            ``'error'``: to get the path corresponding to the error bars.

        only_able : `bool`
            Plot only the contact points that are able to be used in the fit.
            If ``segment='error'`` it will show only the contact points able. If
            segment is any other, the path will be plotted only if both
            immersion and emersion are able, or it is a negative chord.

        ax : `matplotlib.pyplot.Axes`
            The axes where to make the plot. If None, it will use the default axes.

        linestyle : `str`
            Default linestyle used in `matplotlib.pyplot.plot`.
            The difference is that now it accepts ``linestyle='exposure'``, where
            the plot will be a dashed line corresponding to each exposure.
            The blank space between the lines can be interpreted as 'dead time'.

        **kwargs
            Any other kwarg will be parsed directly by `maplotlip.pyplot.plot`.
            The only difference is that the default linewidth ``lw=2``.

        """
        import matplotlib.pyplot as plt

        ax = ax or plt.gca()
        ax.set_xlabel('f (km)')
        ax.set_ylabel('g (km)')
        ax.axis('equal')
        if 'chordlist' not in self._shared_with:
            raise ValueError('{} must be associated to an Occultation to use this function'.format(self.__class__.__name__))

        exposure = False
        if linestyle == 'exposure':
            exposure = True
            kwargs['linestyle'] = '-'
        else:
            kwargs['linestyle'] = linestyle

        kwargs['lw'] = kwargs.get('lw', 2)

        steps = {True: 'exposure', False: 1}
        vals = self.path(segment=segment, step=steps[exposure])
        label = kwargs.pop('label', None)
        var = []
        if only_able:
            if segment == 'error':
                if self.is_able['immersion']:
                    immersion = self.lightcurve.immersion
                    immersion_err = self.lightcurve.immersion_err * u.s
                    var += ax.plot(*self.get_fg(time=[immersion - immersion_err, immersion + immersion_err]), **kwargs)
                if self.is_able['emersion']:
                    emersion = self.lightcurve.emersion
                    emersion_err = self.lightcurve.emersion_err * u.s
                    var += ax.plot(*self.get_fg(time=[emersion - emersion_err, emersion + emersion_err]), **kwargs)
            else:
                if self.is_able.get('immersion') in [None, True] and self.is_able.get('emersion') in [None, True]:
                    var += ax.plot(*vals, **kwargs)
        else:
            var += ax.plot(*vals, **kwargs)
        if len(var) > 0:
            var[0].set_label(label)

    def get_impact_param(self, center_f=0, center_g=0, verbose=True):
        """Gets the impact parameter, minimal distance between the chord and the
        centre position.

        This Chord object must be associated to an Occultation to work, since it
        needs the position of the star and an ephemeris.

        Parameters
        ----------
        center_f : `int`, `float`, default=0
            The coordinate in f of the ellipse center;

        center_g : `int`, `float`, default=0
            The coordinate in g of the ellipse center.

        verbose : `bool`, default=True
            If True, prints the obtained values.


        Returns
        -------
        impact, sense : `list`
            The impact parameter (in km) and the direction of the chord relative
            the ellipse center, North (N), South (S), East (E) and West (W).
        """
        f, g = self.path(segment='full')
        r = np.sqrt((f - center_f)**2 + (g - center_g)**2)
        impact = r.min()
        sense = 'NE'
        if g[np.argmin(r)] < center_g:
            sense = sense.replace('N', 'S')
        if f[np.argmin(r)] < center_f:
            sense = sense.replace('E', 'W')
        if verbose:
            print(self.name)
            print('Impact parameter', np.round(impact, 1), sense)
        return impact, sense

    def get_theoretical_times(self, equatorial_radius, center_f=0, center_g=0, oblateness=0, position_angle=0, sigma=0,
                              step=1, verbose=True):
        """Gets the theoretical times and chord size for a given ellipse.

        This Chord object must be associated to an Occultation to work, since
        it needs the position of the star and an ephemeris.

        Parameters
        ----------
        equatorial_radius : `int`, `float`
            The Equatorial radius (semi-major axis) of the ellipse.

        center_f : `int`, `float`, default=0
            The coordinate in f of the ellipse center

        center_g : `int`, `float`, default=0
            The coordinate in g of the ellipse center

        oblateness : `int`, `float`, default=0
            The oblateness of the ellipse.

        position_angle : `int`, `float`, default=0
            The pole position angle of the ellipse in degrees.
            Zero is in the North direction ('g-positive'). Positive clockwise.

        sigma : `int`, `float`
            Uncertainty of the expected ellipse, in km.

        step : `int`, `float`
            Time resolution of the chord, in seconds.

        verbose : `bool`, default=True
            If True, prints the obtained values.


        Returns
        -------
        theory_immersion_time, theory_emersion_time, theory_chord_size : `list`
            The expected immersion time for the given ellipse, the expected
            emersion time for the given ellipse, and the expected chord size
            for the given ellipse.
        """
        from sora.extra import get_ellipse_points

        time_all = Time(np.arange(self.lightcurve.initial_time.jd, self.lightcurve.end_time.jd, step*u.s.to('d')), format='jd')

        f_all, g_all = self.get_fg(time=time_all)
        df_all = (f_all - center_f)
        dg_all = (g_all - center_g)

        r_all = np.sqrt(df_all**2 + dg_all**2)
        cut = r_all < 1.5*equatorial_radius
        df_path = df_all[cut]
        dg_path = dg_all[cut]
        r_path = r_all[cut]
        time = time_all[cut]

        theta_path = np.arctan2(dg_path, df_path)

        f_ellipse, g_ellipse, r_ellipse, theta = get_ellipse_points(theta_path,
                                                                    equatorial_radius=equatorial_radius,
                                                                    center_f=center_f,
                                                                    center_g=center_g,
                                                                    oblateness=oblateness,
                                                                    position_angle=position_angle)

        ev = r_path < r_ellipse + sigma
        if not np.any(ev):
            if verbose:
                print(self.name)
                print('Negative chord \n')
        try:
            imm = time[ev].jd.argmin()
            eme = time[ev].jd.argmax()
            df_chord = [df_path[ev].min(), df_path[ev].max()]
            dg_chord = [dg_path[ev].min(), dg_path[ev].max()]
            theory_chord_size = np.sqrt((df_chord[1]-df_chord[0])**2 + (dg_chord[1]-dg_chord[0])**2)
            theory_immersion_time = time[ev][imm]
            theory_emersion_time = time[ev][eme]
            if verbose:
                print(self.name)
                print('IMMERSION TIME: {} UTC'.format(theory_immersion_time.iso))
                print('EMERSION  TIME: {} UTC'.format(theory_emersion_time.iso))
                print('CHORD SIZE    : {:.3f} km \n'.format(theory_chord_size))
        except:
            theory_chord_size = 0
            theory_immersion_time = None
            theory_emersion_time = None
        return theory_immersion_time, theory_emersion_time, theory_chord_size

    def __repr__(self):
        """String representation of the Chord Class
        """
        return '<{}: {}>'.format(self.__class__.__name__, self.name)

    def __str__(self):
        """String of the Chord Class used in str(obj) or print(obj)
        """
        from sora.observer import Observer

        string = ['-' * 79, self.observer.__str__()]
        if 'chordlist' in self._shared_with and isinstance(self.observer, Observer):
            star = self._shared_with['chordlist']['star']
            occtime = self._shared_with['chordlist']['time']
            coord = star.get_position(time=occtime, observer=self.observer)
            ephem_altaz = self.observer.altaz(self.lightcurve.time_mean, coord)
            string.append('Target altitude: {:.1f} deg'.format(ephem_altaz[0]))
            string.append('Target azimuth:  {:.1f} deg'.format(ephem_altaz[1]))
        string.append('')
        string.append(self.lightcurve.__str__())
        return '\n'.join(string)
