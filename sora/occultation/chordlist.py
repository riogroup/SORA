from .chord import Chord
from sora.config.list import List
import numpy as np
from astropy.table import Table, vstack
import warnings


__all__ = ['ChordList']


class ChordList(List):
    _allowed_types = (Chord,)  # (Chord, AstrometricChord)
    _set_func = 'add_chord'

    def __init__(self, *, star, body, time):
        """Defines a collection of Chord objects associated to an Occultation.

        This object is not supposed to be defined by the user. It will be automaticaly defined in Occultation.

        Parameters:
            star (Star): The Star occulted.
            body (Body): The occulting Body.
            time (Time): The occultation time.
        """
        self._star = star
        self._body = body
        self._time = time
        self._shared_with = {"chord": {"star": self._star.geocentric(time), "ephem": self._body.ephem}}
        self._shared_with['occultation'] = {}

    def add_chord(self, *, name=None, chord=None, observer=None, lightcurve=None):
        """Add a chord to the occultation chord list

        Parameters:
            name (str): The name of the Chord. It must be unique within the list of chords. If not given, it will use
                the name in chord if chord is directly given or the name in the observer object. If the name in Chord
                already exists in the list, the name parameter can be given to update the chord name. The "name"
                parameter is required if given together with observer and lightcurve. It can not be an empty string.
            chord (Chord): A chord instance defined by the user.
            observer (Observer): An Observer instance defined by the user. It must be given together with lightcurve.
            lightcurve (LightCurve): A lightcurve defined by the user. It must be given together with observer.

        Ways to add a chord:
            obj.add_chord(chord) # where the name in chord will be used.
            obj.add_chord(name, chord) # where the name in chord will be replaced by "name".
            obj.add_chord(observer, lightcurve) # where the name in observer will be used.
            obj.add_chord(name, observer, lightcurve)
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
        lightcurve.set_vel(self._shared_with['occultation']["vel"])
        lightcurve.set_dist(self._shared_with['occultation']["dist"])
        lightcurve.set_star_diam(self._shared_with['occultation']["star_diam"])
        try:
            lightcurve.calc_magnitude_drop(mag_star=self._star.mag['G'], mag_obj=self._body.apparent_magnitude(self._time))
        except:
            lightcurve.bottom_flux = 0.0
            warnings.warn('Magnitude drop was not calculated. Using bottom flux as 0.0.')
        return chord

    def remove_chord(self, *, name):
        """Remove a chord from the chord list and disassociate it from the Occultation.

        Parameters:
            name (str): the name of the chord.
        """
        del(self[name]._shared_with['chordlist'])
        del(self[name])

    def enable(self, *, chord=None, time=None):
        """Enable a contact point of the curve to be used in the fit.

        Parameters:
            chord (str): Name of the chord to enable. If None, it applies to all chords.
            time (None, str): if None, it will enable all contact points.
                if 'immersion' or 'emersion', it will enable respective contact point.
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

        Parameters:
            chord (str): Name of the chord to disable. If None, it applies to all chords.
            time (None, str): if None, it will disable all contact points.
                if 'immersion' or 'emersion', it will disable respective contact point.
        """
        n = 0
        for key in self.keys():
            if chord is None or chord == key:
                self[key].disable(time=time)
                n += 1
        if n == 0:
            raise ValueError('No Chord with name {} is available.'.format(chord))

    def plot_chords(self, *, segment='standard', ignore=None, ax=None, linestyle='-', **kwargs):
        """Plots the on-sky path of this chord.

        Parameters:
            segment (str): The segment to plot the chord. The available options are:
                - 'positive' to get the path between the immersion and emersion times if the chord is positive.
                - 'negative' to get the path between the start and end of observation if the chord is negative.
                - 'standard' to get the 'positive' path if the chord is positive or 'negative' if the chord is negative.
                - 'full' to get  the path between the start and end of observation independent
                    if the chord is positive or negative.
                - 'outer' to get the path outside the 'positive' path, for instance between the start and immersion times
                    and between the emersion and end times.
                - 'error' to get the path corresponding to the error bars.
            ignore (str, list): Name of chord or list of names to ignore in the plot.
            ax (matplotlib.Axes): The axes where to make the plot. If None, it will use the default axes.
            linestyle (str): Default linestyle used in matplotlib.pyplot.plot. The difference is that now it accept
                linestyle='exposure', where the plot will be a dashed line corresponding to each exposure. The blank
                space between the lines can be interpreted as 'dead time'.
            kwargs: Any other kwarg will be parsed directly by maplotlip.pyplot.plot. The only difference is that
                the default linewidth lw=2.

        Returns:
            Default list of plots made by matplotlib.
        """
        n = 0
        keys = list(self.keys())
        plots = []
        if ignore is not None:
            ignore = np.array(ignore, ndmin=1)
        for i in range(len(self)):
            if ignore is not None and keys[i] in ignore:
                continue
            if segment != 'error':
                kwargs['label'] = keys[i]
            try:
                p = self[i].plot_chord(segment=segment, ax=ax, linestyle=linestyle, **kwargs)
                plots += p
            except ValueError:
                n += 1
        if n == len(self):
            warnings.warn('Segment "{}" was not found on any chord'.format(segment))
        return plots

    def summary(self):
        """Prints a table with the summary of the chords.
        """
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
