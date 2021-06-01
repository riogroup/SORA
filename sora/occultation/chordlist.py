import warnings

import numpy as np

from sora.config.list import List
from .chord import Chord

__all__ = ['ChordList']


class ChordList(List):
    """Defines a collection of Chord objects associated to an Occultation.

    This object is not supposed to be defined by the user. It will be automatically
    defined in Occultation.

    Attributes
    ----------
    star : `sora.Star`
        The Star occulted.

    body : `sora.Body`
        The occulting Body.

    time : `astropy.time.Time`
        The occultation time.

    """
    _allowed_types = (Chord,)  # (Chord, AstrometricChord)
    _set_func = 'add_chord'

    def __init__(self, *, star, body, time):

        super().__init__()
        self._star = star
        self._body = body
        self._time = time
        self._shared_with = {"chord": {"star": self._star, "ephem": self._body.ephem, "time": self._time},
                             'occultation': {}}
        self._method_value = 'geocenter'

    def add_chord(self, *, name=None, chord=None, observer=None, lightcurve=None):
        """Add a chord to the occultation chord list

        Parameters
        ----------
        name : `str`
            The name of the Chord. It must be unique within the list of chords.
            If not given, it will use the name in chord if chord is directly
            given or the name in the observer object. If the name in Chord already
            exists in the list, the name parameter can be given to update the
            chord name. The `name` parameter is required if given together with
            observer and lightcurve. It can not be an empty string.

        chord `sora.occultation.Chord`
            A chord instance defined by the user.

        observer : `sora.Observer`
            An Observer instance defined by the user. It must be given together
            with a lightcurve.

        lightcurve : `sora.LightCurve`
            A lightcurve defined by the user. It must be given together with observer.

        Examples
        --------
        Ways to add a chord:

        >>> obj.add_chord(chord) # where the name in chord will be used.

        >>> obj.add_chord(name, chord) # where the name in chord will be replaced by "name".

        >>> obj.add_chord(observer, lightcurve) # where the name in observer will be used.

        >>> obj.add_chord(name, observer, lightcurve)
        """
        if chord and (observer or lightcurve):
            raise ValueError("User must give only chord or (name and observer and lightcurve)")
        if chord:
            name = name or chord.name
        elif observer and lightcurve:
            for key in self.keys():
                if self[key].lightcurve is lightcurve:
                    raise ValueError('lightcurve is already associated to the chord {}'.format(key))
            name = name or observer.name
            chord = Chord(name=name, observer=observer, lightcurve=lightcurve)
        else:
            raise ValueError("User must give chord or (name and observer and lightcurve)")
        self._add_item(name=name, item=chord)
        chord._name = name
        chord._shared_with["chordlist"] = self._shared_with["chord"]
        chord.lightcurve.set_vel(self._shared_with['occultation']["vel"])
        chord.lightcurve.set_dist(self._shared_with['occultation']["dist"])
        chord.lightcurve.set_star_diam(self._shared_with['occultation']["star_diam"])
        try:
            chord.lightcurve.calc_magnitude_drop(mag_star=self._star.mag['G'], mag_obj=self._body.apparent_magnitude(self._time))
        except:
            chord.lightcurve.bottom_flux = 0.0
            warnings.warn('Magnitude drop was not calculated. Using bottom flux as 0.0.')
        chord._method = self._method
        return chord

    # This attribute is to modify the way to calculate f and g on chord.
    @property
    def _method(self):
        return self._method_value

    @_method.setter
    def _method(self, value):
        if value not in ['geocenter', 'observer']:
            raise ValueError('method must be "geocenter" or "observer"')
        self._method_value = value
        for name, chord in self.items():
            chord._method = value

    def remove_chord(self, *, name):
        """Remove a chord from the chord list and disassociate it from the Occultation.

        Parameters
        ----------
        name : `str`
            The name of the chord.
        """
        del(self[name]._shared_with['chordlist'])
        del(self[name])

    @property
    def is_able(self):
        return {name: chord.is_able for name, chord in self.items()}

    def enable(self, *, chord=None, time=None):
        """Enable a contact point of the curve to be used in the fit.

        Parameters
        ----------
        chord : `str`
            Name of the chord to enable. If ``chord=None``, it applies to all chords.

        time : Non`, `str`
            If ``time=None``, it will enable all contact points.
            If 'immersion' or 'emersion', it will enable respective contact point.
        """
        n = 0
        for key in self.keys():
            if chord is None or chord == key:
                self[key].enable(time=time)
                n += 1
        if n == 0:
            raise ValueError('No Chord with name {} is available.'.format(chord))

    def disable(self, *, chord=None, time=None):
        """Disable a contact point of the curve to be used in the fit.

        Parameters
        ----------
        chord : `str`
            Name of the chord to disable. If ``chord=None``, it applies to all chords.

        time : None, `str`
            If ``time=None``, it will disable all contact points.
            If 'immersion' or 'emersion', it will disable respective contact point.
        """
        n = 0
        for key in self.keys():
            if chord is None or chord == key:
                self[key].disable(time=time)
                n += 1
        if n == 0:
            raise ValueError('No Chord with name {} is available.'.format(chord))

    def plot_chords(self, *, segment='standard', ignore_chords=None, only_able=False, ax=None, linestyle='-', **kwargs):
        """Plots the on-sky path of this chord.

        Parameters
        ----------
        segment : `str`
            The segment to plot the chord. The available options are:

            ``'positive'`` to get the path between the immersion and emersion
            times if the chord is positive.

            ``'negative'`` to get the path between the start and end of
            observation if the chord is negative.

            ``'standard'`` to get the 'positive' path if the chord is positive
            or 'negative' if the chord is negative.

            ``'full'`` to get  the path between the start and end of observation
            independent if the chord is positive or negative.

            ``'outer'`` to get the path outside the 'positive' path, for instance
            between the start and immersion times and between the emersion and
            end times.

            ``'error'`` to get the path corresponding to the error bars.

        ignore_chords : `str`, `list`
            Name of chord or list of names to ignore in the plot.

        only_able : `bool`
            Plot only the chords or contact points that are able to be used in the fit.
            If ``segment='error'`` it will show only the contact points able.
            If segment is any other, the path will be plotted only if both
            immersion and emersion are able, or it is a negative chord.

        ax : `matplotlib.pyplot.Axes`
            The axes where to make the plot. If None, it will use the default axes.

        linestyle : `str`
            Default linestyle used in `matplotlib.pyplot.plot`. The difference
            is that now it accept ``linestyle='exposure'``, where the plot will
            be a dashed line corresponding to each exposure. The blank space
            between the lines can be interpreted as 'dead time'.

        **kwargs
            Any other kwarg will be parsed directly by `maplotlip.pyplot.plot`.
            The only difference is that the default linewidth ``lw=2``.

        """
        n = 0
        keys = list(self.keys())
        if ignore_chords is not None:
            ignore_chords = np.array(ignore_chords, ndmin=1)
        for i in range(len(self)):
            if ignore_chords is not None and keys[i] in ignore_chords:
                continue
            if segment != 'error':
                kwargs['label'] = keys[i]
            try:
                _ = self[i].plot_chord(segment=segment, only_able=only_able, ax=ax, linestyle=linestyle, **kwargs)
            except ValueError:
                n += 1
        if n == len(self):
            warnings.warn('Segment "{}" was not found on any chord'.format(segment))

    def summary(self):
        """Prints a table with the summary of the chords.
        """
        from astropy.table import Table, vstack

        tables = []
        for key in self.keys():
            tt = Table()
            obs = self[key].observer
            lc = self[key].lightcurve
            max_val = 0
            cols = []
            colnames = ['Name', 'Longitude', 'Latitude', 'status', 'time', 'f', 'g']
            itens = [key, obs.lon.to_string(), obs.lat.to_string()]
            row = []
            times = {'Initial Time': 'initial_time', 'Immersion': 'immersion', 'Emersion': 'emersion', 'End Time': 'end_time'}
            for i in times:
                val = getattr(lc, times[i], None)
                if val is not None:
                    row.append([i, val.iso, *['{:.2f}'.format(n) for n in self[key].get_fg(time=val)]])
            for t in np.array(row).T:
                itens.append(t.tolist())
            for item in itens:
                v = np.array(item, ndmin=1).tolist()
                cols.append(v)
                if len(v) > max_val:
                    max_val = len(v)
            for i in range(len(cols)):
                if len(cols[i]) < max_val:
                    for n in range(max_val - len(cols[i])):
                        cols[i].append('')
                tt[colnames[i]] = cols[i]
            tables.append(tt)
        tabela = vstack(tables)
        tabela.pprint_all()

    def get_impact_param(self, chords='all_chords', center_f=0, center_g=0, verbose=True):
        """Get the impact parameter, minimal distance between the chord and the
        centre position.

        This Chord object must be associated to an Occultation to work, since it
        needs the position of the star and an ephemeris.

        Parameters
        ----------
        chords : `int`, `str`, default='all_chords'
            Index or names of the chords to be considered.

        center_f : `int`, `float`, default=0
            The coordinate in f of the ellipse center.

        center_g : `int`, `float`, default=0
            The coordinate in g of the ellipse center.

        verbose : `bool`
            If True, prints the obtained values.


        Returns
        -------
        impact, sense, chord_name : `list`
            The impact parameter (in km), the direction of the chord relative
            the ellipse center, North (N), South (S), East (E) and West (W), and
            the name of the chord
        """
        impact = np.array([])
        sense = np.array([])
        names = np.array([])
        if chords == 'all_chords':
            chords = range(len(self))
        for i in chords:
            chord = self[i]
            im, se = chord.get_impact_param(center_f=center_f, center_g=center_g, verbose=verbose)
            impact = np.append(impact, im)
            sense = np.append(sense, se)
            names = np.append(names, chord.name)
        return impact, sense, names

    def get_theoretical_times(self, equatorial_radius, chords='all_chords', center_f=0, center_g=0, oblateness=0,
                              position_angle=0, sigma=0, step=1, verbose=True):
        """Get the theoretical times and chords sizes for a given ellipse.

        This Chord object must be associated to an Occultation to work, since it needs
        the position of the star and an ephemeris.

        Parameters
        ----------
        chords : `int`, `str`, default='all_chords'
            Index or names of the chords to be considered.

        equatorial_radius : `int`, `float`
            The Equatorial radius (semi-major axis) of the ellipse.

        center_f : `int`, `float`, default=0
            The coordinate in f of the ellipse center.

        center_g : `int`, `float`, default=0
            The coordinate in g of the ellipse center.

        oblateness : `int`, `float`, default=0
            The oblateness of the ellipse.

        position_angle : `int`, `float`, default=0
            The pole position angle of the ellipse in degrees.
            Zero is in the North direction ('g-positive'). Positive clockwise.

        sigma : `int`, `float`
            Uncertainty of the expected ellipse, in km.

        step : `int`, `float`
            Time resolution of the chord, in seconds.

        verbose : `bool`
            If True, prints the obtained values.


        Returns
        -------
        theory_immersion_time, theory_emersion_time, theory_chord_size, chord_name : `list`
            The expected immersion time for the given ellipse, the expected
            emersion time for the given ellipse, the expected chord size for the
            given ellipse, and the name of the chord.
        """
        theory_immersion_time = np.array([])
        theory_emersion_time = np.array([])
        theory_chord_size = np.array([])
        names = np.array([])
        if chords == 'all_chords':
            chords = range(len(self))
        for i in chords:
            chord = self[i]
            tit, tet, tcs = chord.get_theoretical_times(equatorial_radius=equatorial_radius, center_f=center_f,
                                                        center_g=center_g, oblateness=oblateness, position_angle=position_angle,
                                                        sigma=sigma, step=step, verbose=verbose)
            theory_immersion_time = np.append(theory_immersion_time, tit)
            theory_emersion_time = np.append(theory_emersion_time, tet)
            theory_chord_size = np.append(theory_chord_size, tcs)
            names = np.append(names, chord.name)
        return theory_immersion_time, theory_emersion_time, theory_chord_size, names
