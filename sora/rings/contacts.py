#sora/rings/contacts.py
from __future__ import annotations
from dataclasses import dataclass
from typing import Any, Dict
import matplotlib.pyplot as plt
import astropy.units as u
import numpy as np


@dataclass
class RingContact:
    ring : object
    label: str
    chi2: Any
    chord: Any
    contact: str

    time_mean: float | None

    immersion: float | None
    immersion_err: float | None

    emersion: float | None
    emersion_err: float | None

    opacity: float | None
    opacity_err: float | None

    properties: Dict[str, Any] | None = None

    def plot_contact(self, *, segment='central_point', plot_ring=False, 
                      ring_radius=None, ax=None, time_direction=False, 
                      **kwargs):    
    
        ax = ax or plt.gca()
        ax.set_xlabel('f (km)')
        ax.set_ylabel('g (km)')
        ax.axis('equal')
        
        dur = (self.emersion - self.immersion).sec*u.s        
            
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
            ax.plot(*self.chord.get_fg(time=[self.time_mean - dur/2 - self.immersion_err*u.s,
                                                    self.time_mean + dur/2 + self.emersion_err*u.s]),
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

class RingContactList(dict):
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

    # ------------------------------------------------------------
    def __getitem__(self, key: str | int):
        """
        Access by label (string) or by integer index.
        """
        if isinstance(key, str):
            return super().__getitem__(key)
        elif isinstance(key, int):
            return list(self.values())[key]
        else:
            raise TypeError(
                "Key must be a string (label) or an integer index."
            )

    # ------------------------------------------------------------
    def __repr__(self):
        """
        Representation, similar to ChordList.
        """
        out = ["<RingOccultationList:"]
        for i, occ in enumerate(self.values()):
            out.append(f"    {i}: RingOccultation('{occ.label}')")
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
