"""
Unified Detection Limits module for SORA
---------------------------------------

This file replaces the old lightcurve.detectionlimits module and also
handles chord-based projected/diffracted limits. The same class is used for:

    - LightCurve  → apparent (sky-plane) detection limits
    - Chord       → projected (ring-plane corrected) limits

The apparent noise model is exactly the same as the implementation provided
previously for LightCurve.

For Chord, the apparent limits are computed from the associated LightCurve and
then optionally converted into physical/projected limits when an orientation
is supplied (P, B, pole, or ring).
"""

import numpy as np
import warnings
import astropy.units as u
from astropy.stats import sigma_clip

warnings.simplefilter("always", UserWarning)

__all__ = ["DetectionLimits"]


class DetectionLimits:
    """
    Unified detection limits for LightCurve and Chord.

    Parameters
    ----------
    source : LightCurve or Chord
        Object from which the detection limits are derived.
    P : astropy.units.Quantity, optional
        Position angle of the ring plane (degrees), measured east of north.
        Used only when `source` is a Chord.
    B : astropy.units.Quantity, optional
        Opening angle of the ring plane (degrees); B = 0° is edge-on.
        Used only when `source` is a Chord.
    pole : astropy.coordinates.SkyCoord, optional
        Pole orientation of the ring. Used only for information in the output
        string (it is not used to compute P and B inside this class).
    ring : sora.rings.core.Ring, optional
        Ring object. If given and P/B are not provided, this class will call
        ``ring.get_ring_orientation(time=occtime)`` using the chord's closest
        approach time stored in ``chord._shared_with['chordlist']['time']``.

    Notes
    -----
    LightCurve:
        - Only apparent opacity, optical depth and equivalent width are used.
        - No geometry or projection is applied.

    Chord:
        - The same apparent limits are computed from the underlying LightCurve.
        - If no orientation (P,B,pole,ring) is provided, only apparent limits
          are available; projected quantities (normal opacity, etc.) will be
          returned as None, and the string output will indicate that no pole/
          orientation is available.
        - If orientation is provided, Fresnel correction and projection via
          sin(B) are applied.
    """

    def __init__(self, source, *, P=None, B=None, pole=None, ring=None):
        self.source = source
        self.kind = self._infer_kind(source)

        # Shared attributes
        self.flux = None
        self.flux_err = None

        # Geometry / orientation (Chord only)
        self._P = None
        self._B = None
        self._used_pole = None

        # Initialize per type
        if self.kind == "lightcurve":
            self._init_from_lightcurve()

        elif self.kind == "chord":
            self._init_from_chord(P=P, B=B, pole=pole, ring=ring)

        else:
            raise TypeError(f"DetectionLimits does not support type {type(source)}")


    def _infer_kind(self, src):
        from sora.lightcurve.core import LightCurve
        from sora.occultation.chord import Chord

        if isinstance(src, LightCurve):
            return "lightcurve"
        if isinstance(src, Chord):
            return "chord"
        return None


    @staticmethod
    def _compute_flux_from_lightcurve(lc):
        """
        Prepare flux and flux_err arrays for detection limits, masking the
        occultation interval (±1 exposure) if immersion/emersion are defined.
        """
        flux = getattr(lc, "flux", None)
        flux_err = getattr(lc, "dflux", None)

        if flux is None or lc.time is None:
            return None, None

        immersion = getattr(lc, "immersion", None)
        emersion = getattr(lc, "emersion", None)
        tref = getattr(lc, "tref", None)

        flux = np.array(flux)

        if immersion and emersion and tref:
            imm = (immersion - tref).sec - lc.exptime
            eme = (emersion - tref).sec + lc.exptime
            mask = (lc.time < imm) + (lc.time > eme)
            flux = flux[mask]

            if flux_err is not None:
                flux_err = np.array(flux_err)[mask]
        else:
            if flux_err is not None:
                flux_err = np.array(flux_err)

        if flux_err is not None and len(flux_err) != len(flux):
            flux_err = None

        return flux, flux_err


    def _init_from_lightcurve(self):
        """
        Sets up the inputs for the apparent detection-limit code.
        """
        lc = self.source
        self.flux, self.flux_err = self._compute_flux_from_lightcurve(lc)

        # No geometry in LightCurve mode
        self._P = None
        self._B = None
        self._used_pole = None


    def _init_from_chord(self, *, P=None, B=None, pole=None, ring=None):
        """
        In chord mode, the apparent limits come from the underlying LightCurve,
        and then we optionally apply:

            - Fresnel correction (Cuzzi 1985)
            - projection using sin(B), if B is provided
        """
        chord = self.source
        lc = chord.lightcurve

        # Apparent part: reuse the LC logic
        self.flux, self.flux_err = self._compute_flux_from_lightcurve(lc)

        used_pole = None
        occtime = None

        # Get occultation time for ring orientation
        if hasattr(chord, "_shared_with"):
            occtime = chord._shared_with.get("chordlist", {}).get("time", None)

        # CASE 1 — ring provided → get P,B from ring
        if ring is not None and occtime is not None and (P is None or B is None):
            try:
                P_ring, B_ring = ring.get_ring_orientation(time=occtime)
                if P is None:
                    P = P_ring
                if B is None:
                    B = B_ring

                # Try to build a pole string from ring geometry if available
                if hasattr(ring, "geometry") and hasattr(ring.geometry, "pole"):
                    try:
                        used_pole = ring.geometry.pole.to_string("hmsdms")
                    except Exception:
                        used_pole = "Ring geometry pole"
                elif hasattr(ring, "pole_orientation"):
                    try:
                        used_pole = ring.pole_orientation.icrs.to_string("hmsdms")
                    except Exception:
                        used_pole = "Ring pole_orientation"
            except Exception:
                # fall back silently; user may have provided P,B directly
                pass

        # CASE 2 — pole provided only for info (we do NOT compute P,B from it here)
        if pole is not None:
            try:
                used_pole = pole.icrs.to_string("hmsdms")
            except Exception:
                used_pole = str(pole)

        # If user provided P,B explicitly and no pole
        if (P is not None or B is not None) and used_pole is None:
            used_pole = "User orientation"

        self._P = P
        self._B = B
        self._used_pole = used_pole


    def _estimate_noise(self):
        """Same algorithm used originally in LightCurve detectionlimits."""
        if self.flux is None or len(self.flux) == 0:
            return None

        clipped = sigma_clip(self.flux, sigma=5, maxiters=5)
        clean_flux = clipped.data[~clipped.mask]

        if len(clean_flux) < 10:
            clean_flux = self.flux

        if self.flux_err is not None:
            clean_err = (
                self.flux_err[~clipped.mask]
                if len(self.flux_err) == len(self.flux)
                else self.flux_err
            )
            weights = 1 / clean_err**2
            wmean = np.average(clean_flux, weights=weights)
            wvar = np.average((clean_flux - wmean) ** 2, weights=weights)
            return np.sqrt(wvar)

        return np.std(clean_flux, ddof=1)

    def apparent_opacity(self, sigma=1):
        """
        Apparent opacity detection limit (sky-plane), same for LightCurve and Chord.
        """
        noise = self._estimate_noise()
        if noise is None:
            return None
        return noise * sigma

    def apparent_optical_depth(self, sigma=1):
        """
        Apparent optical depth detection limit: τ = -ln(1 - sigma * noise).
        """
        noise = self._estimate_noise()
        if noise is None:
            return None
        return -np.log(1 - sigma * noise)

    def apparent_equivalent_width(self, sigma=1):
        """
        Apparent equivalent width detection limit (in km).
        """
        # Determine the underlying LightCurve
        if self.kind == "lightcurve":
            lc = self.source
        else:
            lc = self.source.lightcurve  # chord mode

        noise = self._estimate_noise()
        if noise is None:
            return None

        return sigma * noise * lc.vel * lc.exptime

    # ------------------------------------------------------------------
    # PHYSICAL / PROJECTED LIMITS (CHORD ONLY)
    # ------------------------------------------------------------------
    def opacity(self, sigma=1):
        """
        Physical opacity corrected for Fresnel diffraction (Chord only).

        For LightCurve mode, this simply returns the apparent opacity.
        """
        # LightCurve → return apparent only
        if self.kind == "lightcurve":
            return self.apparent_opacity(sigma=sigma)

        # Chord
        opa_app = self.apparent_opacity(sigma=sigma)
        if opa_app is None:
            return None

        lc = self.source.lightcurve

        # Fresnel correction following Cuzzi (1985), Roques et al. (1987)
        particle_size = 0.001  # km (assumed)
        wavelength = lc.lambda_0 * u.micrometer.to("km")
        dist_km = lc.dist * u.au.to("km")
        spatial_resolution = abs(lc.exptime * lc.vel)
        airy_scale = (wavelength / (2 * particle_size)) * dist_km

        if airy_scale > spatial_resolution:
            return 1 - np.sqrt(1 - opa_app)
        else:
            return opa_app

    def optical_depth(self, sigma=1):
        """
        Optical depth τ = -ln(1 - opacity).
        """
        op = self.opacity(sigma=sigma)
        if op is None:
            return None
        return -np.log(1 - op)

    def normal_opacity(self, sigma=1):
        """
        Normal opacity, projected using sin(|B|).

        Returns None if B is not available.
        """
        if self._B is None:
            return None
        op = self.opacity(sigma=sigma)
        if op is None:
            return None
        return op * abs(np.sin(self._B).value)

    def normal_optical_depth(self, sigma=1):
        """
        Normal optical depth, projected using sin(|B|).

        Returns None if B is not available.
        """
        if self._B is None:
            return None
        tau = self.optical_depth(sigma=sigma)
        if tau is None:
            return None
        return tau * abs(np.sin(self._B).value)


    def __str__(self):
        # LightCurve: only apparent limits
        if self.kind == "lightcurve":
            if self.flux is None or len(self.flux) == 0:
                return "\nNo flux data available to compute detection limits.\n"

            return (
                "\nDetection limits (3-sigma):\n"
                f"    Apparent opacity:             {self.apparent_opacity(sigma=3):.3f}\n"
                f"    Apparent optical depth:       {self.apparent_optical_depth(sigma=3):.3f}\n"
                f"    Apparent equivalent width:    {self.apparent_equivalent_width(sigma=3):.3f} km\n\n"
            )

        # Chord: check if we have orientation
        if self.kind == "chord":
            if self._B is None:
                return (
                    "\nPole/orientation not available to compute projected detection limits.\n"
                    "Use chord.detection_limits(P=..., B=...) or chord.detection_limits(ring=...) "
                    "to obtain projected quantities.\n"
                )

            return (
                "\nProjected detection limits (3-sigma):\n"
                f"    Reference pole:               {self._used_pole}\n"
                f"    Opening angle (B):            {self._B:.3f}\n"
                f"    Position angle (P):           {self._P:.3f}\n"
                f"    Normal opacity:               {self.normal_opacity(sigma=3):.3f}\n"
                f"    Normal optical depth:         {self.normal_optical_depth(sigma=3):.3f}\n"
            )

        return ""
