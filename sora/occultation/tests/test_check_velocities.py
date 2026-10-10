import re
from types import SimpleNamespace

import numpy as np
import pytest

from sora.extra.utils import get_ellipse_points
from sora.occultation import Occultation


@pytest.mark.parametrize('position_angle', [0.0, 60.0])
def test_check_velocities_uses_position_angle(position_angle, capsys):
    params = {'center_f': [100.0, 0.0], 'center_g': [-50.0, 0.0], 'equatorial_radius': [800.0, 0.0],
              'oblateness': [0.4, 0.0], 'position_angle': [position_angle, 0.0]}
    kw = dict(equatorial_radius=800.0, oblateness=0.4, center_f=100.0, center_g=-50.0,
              position_angle=position_angle)
    velocity = np.array([20.0, -5.0])
    vals, expected = {}, []
    for key, theta in ((1, 0.6), (2, 3.7)):
        eps = 1e-6  # normal from finite differences of the ellipse, independent of the code under test
        (x0, y0), (xp, yp), (xm, ym) = (get_ellipse_points(t, **kw)[:2] for t in (theta, theta + eps, theta - eps))
        tangent = np.array([xp - xm, yp - ym]) / np.hypot(xp - xm, yp - ym)
        vals[key] = np.array([x0, y0, *velocity])
        expected.append(abs(np.dot([-tangent[1], tangent[0]], velocity)))
    chord = SimpleNamespace(lightcurve=SimpleNamespace(vel=np.hypot(*velocity), immersion=1, emersion=2),
                            get_fg=lambda time=None, vel=False: vals[time])

    Occultation.check_velocities(SimpleNamespace(fitted_params=params, chords={'chord': chord}))

    printed = [float(v) for v in re.findall(r'Radial Velocity: ([0-9.]+)', capsys.readouterr().out)]
    np.testing.assert_allclose(printed, expected, atol=2e-3)
