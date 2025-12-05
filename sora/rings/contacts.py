#sora/rings/contacts.py
from __future__ import annotations
from dataclasses import dataclass
from typing import Any, Dict
import matplotlib.pyplot as plt
import astropy.units as u
import numpy as np
from collections import OrderedDict


class RingContact:
    """
    Representa um contato de anel detectado em uma curva de luz.
    Pode ser inicializado a partir de:
        - um modelo SORA
        - um objeto de chi2
        - valores numéricos # TODO ou Time.
    """

    def __init__(self, *, ring, label, chord, contact,
                 sigma=1, model=None, chi2=None,
                 immersion=None, emersion=None,
                 immersion_err=None, emersion_err=None):

        self.ring = ring
        self.label = label
        self.chord = chord
        self.contact = contact

        self.model = model
        self.chi2 = chi2
        self._sigma = sigma
        self._immersion_num = immersion
        self._emersion_num = emersion
        self._immersion_err = immersion_err
        self._emersion_err = emersion_err

        modes = [model is not None,
                 chi2 is not None,
                 immersion is not None or emersion is not None]
        if sum(modes) == 0:
            raise ValueError("Must provide model OR chi2 OR numeric values")
        if sum(modes) > 1:
            raise ValueError("Provide only ONE type of input")

    @property
    def immersion(self):
        if self._immersion_num is not None:
            return self.chord.lightcurve.tref + self._immersion_num*u.s
        if self.model is not None:
            return self.chord.lightcurve.tref + self.model.params.get("immersion")*u.s
        if self.chi2 is not None:
            return self.chord.lightcurve.tref + self.chi2.get_nsigma(sigma=self._sigma)['immersion'][0]*u.s
        return None

    @property
    def emersion(self):
        if self._emersion_num is not None:
            return self._emersion_num
        if self.model is not None:
            return self.chord.lightcurve.tref + self.model.params.get("emersion")*u.s

        if self.chi2 is not None:
            return self.chord.lightcurve.tref + self.chi2.get_nsigma(sigma=self._sigma)['emersion'][0]*u.s
        return None
    
    @property
    def immersion_err(self):
        if self._immersion_err is not None:
            return self.immersion_err
        if self.model is not None:
            return float(0.0) # TODO - ajustar o model.params para ter o _err também
        if self.chi2 is not None:
            return self.chi2.get_nsigma(sigma=self._sigma)['immersion'][1]
        return None
    
    @property
    def emersion_err(self):
        if self._immersion_err is not None:
            return self.emersion_err
        if self.model is not None:
            return float(0.0) # TODO - ajustar o model.params para ter o _err também
        if self.chi2 is not None:
            return self.chi2.get_nsigma(sigma=self._sigma)['emersion'][1]
        return None

    @property
    def time_mean(self):
        i = (self.immersion - self.chord.lightcurve.tref).sec
        e = (self.emersion - self.chord.lightcurve.tref).sec
        if i is None or e is None:
            return None
        return self.chord.lightcurve.tref + (0.5*(i+e))*u.s

    def plot_contact(self, *, segment='central_point', plot_ring=False, 
                      ring_radius=None, ax=None, time_direction=False, 
                      **kwargs):    
    
        ax = ax or plt.gca()
        ax.set_xlabel('f (km)')
        ax.set_ylabel('g (km)')
        ax.axis('equal')
                    
        var = [] 
        
        if segment == 'standard':
            ax.plot(*self.chord.get_fg(time=[self.immersion - self.immersion_err*u.s,
                                            self.immersion + self.immersion_err*u.s]),
                    color='r', linewidth=2, zorder=1)
            ax.plot(*self.chord.get_fg(time=[self.emersion - self.emersion_err*u.s,
                                            self.emersion + self.emersion_err*u.s]),
                    color='r', linewidth=2, zorder=1)
            ax.plot(*self.chord.get_fg(time=[self.immersion, self.emersion]),
                    **kwargs, linewidth=1, zorder=5)

        if segment == 'central_point':
            fc, gc = self.chord.get_fg(time=[self.time_mean])
            dur = (self.emersion - self.immersion).sec
            ax.plot(*self.chord.get_fg(time=[self.time_mean - (dur/2 - self.immersion_err)*u.s,
                                                    self.time_mean + (dur/2 + self.emersion_err)*u.s]),
                        color='r', linewidth=2, zorder=1)
            ax.scatter(fc, gc, marker='.', zorder=5)
                        
        if plot_ring:
            if ring_radius is None:
                ring_radius = getattr(self.ring, "radius", None)
                if getattr(ring_radius, "value", None) is not None:
                    ring_radius = ring_radius.value

            if ring_radius is not None:
                from sora.extra import draw_ellipse
                P, B = self.ring.get_ring_orientation(time=self.chord.lightcurve.tref)
                draw_ellipse(equatorial_radius=ring_radius, 
                            position_angle=P.value, 
                            oblateness=1-abs(np.sin(B)),
                            ax=ax, zorder=0, lw=1
                            )
            else:
                warnings.warn(
                    f"Ring '{self.ring.name}' initialized without radius information. Please furnish `radius`.",
                    UserWarning
                    )


        if len(var) > 0:
                if time_direction is True:
                    if segment != 'error':
                        in_fg = self.chord.get_fg(time=[self.chord.lightcurve.initial_time])[0]
                        out_fg = self.chord.get_fg(time=[self.chord.lightcurve.end_time])[0]
                        if out_fg[0] > in_fg[0]:
                            add_arrow(var[0], direction='right')
                        elif out_fg[0] < in_fg[0]:
                            add_arrow(var[0], direction='left')
                        else:
                            pass
    def get_fg(self, vel=False):
        fc, gc, vfc, vgc = self.chord.get_fg(time=[self.time_mean], vel=True)
        if not vel:
            return fc, gc
        return fc, gc, vfc, vgc

class RingContactList(OrderedDict):
    """
    Container for RingContact objects.

    Behaves like a hybrid of:

    - a dictionary   (access by label:  list['Q1R_c1_1st'])
    - a list         (access by index:  list[0])

    Provides automatic validation to prevent duplicated entries,
    ensuring consistency of ring-contact registrations, in analogy
    with the behavior of Occultation.ChordList.

    Notes
    -----
    - Keys are labels (strings).
    - Values are RingContact instances.
    """

    def validate_add(self, occ: RingContact):
        """
        Raises ValueError if the RingContact conflicts with
        one already stored in the container.

        A conflict occurs when:

        1. The same label is already present.
        2. The same chord and same contact ('1st' or '2nd') is registered.
        """

        if occ.label in self:
            raise ValueError(
                f"Ring occultation '{occ.label}' is already registered."
            )

        new_contact = occ.label.split("_")[-1]

        for existing in self.values():
            existing_contact = existing.label.split("_")[-1]
            if (existing.chord is occ.chord) and (existing_contact == new_contact):
                raise ValueError(
                    f"Contact '{new_contact}' already exists for chord "
                    f"'{occ.chord.name}' (label '{existing.label}')."
                )

    def __setitem__(self, key: str, occ: RingContact):
        self.validate_add(occ)
        super().__setitem__(key, occ)

    def __getitem__(self, key: str | int):
        """
        Access by label (string) or by integer index.
        """
        if isinstance(key, str):
            return super().__getitem__(key)
        elif isinstance(key, int):
            return list(self.values())[key]
        else:
            raise TypeError("Key must be string (label) or integer index.")

    def __iter__(self):
        """Iterate over RingContact objects (like ChordList)."""
        return iter(self.values())
    
    def __repr__(self):
        """
        Representation, similar to ChordList.
        """
        out = ["<RingContactList:"]
        for i, occ in enumerate(self.values()):
            out.append(f"    {i}: RingContact({occ.chord.name}, {occ.contact})")
        out.append(">")
        return "\n".join(out)

    def plot_contacts(self, *, segment='central_point', plot_ring=False, 
                      ring_radius=None, ax=None, time_direction=False, 
                      **kwargs):
        for i, occ in enumerate(self.values()):
            self[i].plot_contact(segment=segment, 
                                 plot_ring=plot_ring, 
                                 ring_radius=ring_radius, ax=ax, 
                                 time_direction=time_direction, 
                                 **kwargs)
