import numpy as np
import astropy.units as u
from numpy.testing import assert_allclose

from einsteinpy.coordinates import SphericalDifferential
from einsteinpy.metric import Schwarzschild


def test_f_vec_bl_schwarzschild():
    M = 6.73317655e26 * u.kg
    sph = SphericalDifferential(
        e0=0. * u.s,
        e1=1e6 * u.m,
        e2=4 * np.pi / 5 * u.rad,
        e3=0. * u.rad,
        u0=0. * u.m / u.s,
        u1=0. * u.rad / u.s,
        u2=2e6 * u.rad / u.s
    )
    f_vec_expected = np.array(
        [
            3.92128321e+03, 0.00000000e+00, 0.00000000e+00, 2.00000000e+06,
            -0.00000000e+00, 1.38196394e+18, -1.90211303e+12, -0.00000000e+00
        ]
    )

    ms = Schwarzschild(coords=sph, M=M)
    state = np.hstack((sph.position(), sph.velocity(ms)))

    f_vec = ms._f_vec(0., state)
    f_vec

    assert isinstance(f_vec, np.ndarray)
    assert_allclose(f_vec_expected, f_vec, rtol=1e-8)
