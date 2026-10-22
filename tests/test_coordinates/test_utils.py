"""
This is a test module, that tests the non-coordinate conversion utilities in ``einsteinpy.coordinates.utils``.

"""
import astropy.units as u
import numpy as np
import pytest
from numpy.testing import assert_allclose

from einsteinpy.coordinates import BoyerLindquistDifferential, SphericalDifferential
from einsteinpy.coordinates.utils import lorentz_factor, v0
from einsteinpy.metric import Schwarzschild, Kerr, KerrNewman

from einsteinpy import constant

_c = constant.c.value


def test_lorentz_factor():
    """
    Tests, if the value of Lorentz Factor, calculated using \
    ``einsteinpy.coordindates.utils.lorentz_factor()`` is the same as, \
    that calculated by brute force

    """
    # Calculated using v0()
    gamma = lorentz_factor(-0.1, -0.01, 0.05)

    # Calculated by brute force
    v_vec = np.array([-0.1, -0.01, 0.05])
    v_norm2 = v_vec[0]**2 + v_vec[1]**2 + v_vec[2]**2
    gamma_b = 1 / np.sqrt(1 - v_norm2 / _c ** 2)

    assert_allclose(gamma, gamma_b, rtol=1e-8)


@pytest.fixture
def sph():
    return SphericalDifferential(
        0. * u.s,
        1. * u.m,
        np.pi / 2 * u.rad,
        0.1 * u.rad,
        0. * u.m / u.s,
        0.1 * u.rad / u.s,
        951 * u.rad / u.s
    )


@pytest.fixture
def bl():
    return BoyerLindquistDifferential(
        0. * u.s,
        1. * u.m,
        np.pi / 2 * u.rad,
        0.1 * u.rad,
        0. * u.m / u.s,
        0.1 * u.rad / u.s,
        951 * u.rad / u.s
    )


def test_v0(sph, bl):
    """
    Tests, if the 4-Velocity in KerrNewman Metric is the same as that in Kerr Metric, \
    in the limit Q -> 0 and if it becomes the same as that in Schwarzschild \
    Metric, in the limits, a -> 0 & Q -> 0

    """
    M = 1e24 * u.kg

    ms = Schwarzschild(coords=sph, M=M)
    mk = Kerr(coords=bl, M=M, a=0.5 * u.one)
    mk0 = Kerr(coords=bl, M=M, a=0. * u.one)
    mkn = KerrNewman(coords=bl, M=M, a=0.5 * u.one, Q=0. * u.C)
    mkn0 = KerrNewman(coords=bl, M=M, a=0. * u.one, Q=0. * u.C)

    v4vec_s = sph.velocity(ms)
    v4vec_k = bl.velocity(mk)
    v4vec_k0 = bl.velocity(mk0)
    v4vec_kn = bl.velocity(mkn)
    v4vec_kn0 = bl.velocity(mkn0)

    assert_allclose(v4vec_s, v4vec_k0, rtol=1e-8)
    assert_allclose(v4vec_k, v4vec_kn, rtol=1e-8)
    assert_allclose(v4vec_kn0, v4vec_s, rtol=1e-8)


def test_compare_vt_schwarzschild():
    """
    Tests, if the value of timelike component of 4-Velocity in Schwarzschild spacetime, \
    calculated using ``einsteinpy.coordindates.utils.v0()`` is the same as that calculated by \
    brute force

    """
    # Calculated using v0()
    M = 1e24 * u.kg
    sph = SphericalDifferential(
        0. * u.s,
        1. * u.m,
        np.pi / 2 * u.rad,
        0.1 * u.rad,
        -0.1 * u.m / u.s,
        -0.01 * u.rad / u.s,
        0.05 * u.rad / u.s
    )

    ms = Schwarzschild(coords=sph, M=M)
    x_vec = sph.position()
    ms_mat = ms.metric_covariant(x_vec)
    v_vec = sph.velocity(ms)

    sph.u0 = (ms,)  # Setting u0
    vt_s = sph.u0  # Getting u0

    # Calculated by brute force
    A = ms_mat[0, 0]
    C = ms_mat[1, 1] * v_vec[1]**2 + ms_mat[2, 2] * v_vec[2]**2 + ms_mat[3, 3] * v_vec[3]**2 - _c ** 2
    D = - 4 * A * C
    vt_sb = np.sqrt(D) / (2 * A)

    assert_allclose(vt_s.value, vt_sb, rtol=1e-8)


@pytest.fixture
def sph2():
    return SphericalDifferential(
        0. * u.s,
        1. * u.m,
        np.pi / 2 * u.rad,
        0.1 * u.rad,
        -0.1 * u.m / u.s,
        -0.01 * u.rad / u.s,
        0.05 * u.rad / u.s
    )


@pytest.fixture
def bl2():
    return BoyerLindquistDifferential(
        0. * u.s,
        1. * u.m,
        np.pi / 2 * u.rad,
        0.1 * u.rad,
        -0.1 * u.m / u.s,
        -0.01 * u.rad / u.s,
        0.05 * u.rad / u.s
    )


def test_compare_vt_schwarzschild_kerr_kerrnewman(sph2, bl2):
    """
    Tests, whether the timelike component of 4-Velocity in KerrNewman Metric is the same as that \
    in Kerr Metric, in the limit Q -> 0 and if it becomes the same as that in Schwarzschild \
    Metric, in the limits, a -> 0 & Q -> 0

    """
    M = 1e24 * u.kg

    ms = Schwarzschild(coords=sph2, M=M)
    mk = Kerr(coords=bl2, M=M, a=0.5 * u.one)
    mk0 = Kerr(coords=bl2, M=M, a=0. * u.one)
    mkn = KerrNewman(coords=bl2, M=M, a=0.5 * u.one, Q=0. * u.C)
    mkn0 = KerrNewman(coords=bl2, M=M, a=0. * u.one, Q=0. * u.C)

    v_vec_ms = sph2.velocity(ms)
    v_vec_mk = bl2.velocity(mk)
    v_vec_mk0 = bl2.velocity(mk0)
    v_vec_mkn = bl2.velocity(mkn)
    v_vec_mkn0 = bl2.velocity(mkn0)

    assert_allclose(v_vec_ms, v_vec_mk0, rtol=1e-8)
    assert_allclose(v_vec_mk, v_vec_mkn, rtol=1e-8)
    assert_allclose(v_vec_mkn0, v_vec_ms, rtol=1e-8)
