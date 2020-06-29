"""
This is a test module, that tests the non-coordinate conversion utilities in ``einsteinpy.coordinates.utils``.

"""

import numpy as np
import pytest
from numpy.testing import assert_allclose

from einsteinpy.metric import Schwarzschild, Kerr, KerrNewman
from einsteinpy.coordinates.utils import lorentz_factor, four_position, four_velocity, v0
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


@pytest.mark.parametrize(
    "t, x_vec",
    [
        (
            0.,
            np.array([306., np.pi / 2, np.pi / 2]),
        ),
        (
            100.,
            np.array([1e3, 0.15, np.pi / 2]),
        ),
    ],
)
def test_four_position(t, x_vec):
    """
    Tests, if the 4-Position, returned by ``einsteinpy.coordindates.utils.four_position()`` \
    is the same as, that calculated manually

    """
    x_4vec = four_position(t, *x_vec)
    x_4vec_b = np.append([_c * t], x_vec)

    assert_allclose(x_4vec, x_4vec_b)


@pytest.mark.parametrize(
    "v_vec, time_like",
    [
        (
            np.array([0., 0.1, 951.]),
            True,
        ),
        (
            np.array([0.2 * _c, 0.5e-5 * _c, 1e-4 * _c]),
            False,
        ),
    ],
)
def test_four_velocity(v_vec, time_like):
    """
    Tests, if the 4-Velocity in KerrNewman Metric is the same as that in Kerr Metric, \
    in the limit Q -> 0 and if it becomes the same as that in Schwarzschild \
    Metric, in the limits, a -> 0 & Q -> 0

    """
    M = 1e24
    x_vec = np.array([1.0, np.pi / 2, 0.1])

    x_vec = np.array([1.0, np.pi / 2, 0.1])

    ms = Schwarzschild(M=M)
    mk = Kerr(coords="BL", M=M, a=0.5)
    mk0 = Kerr(coords="BL", M=M, a=0.)
    mkn = KerrNewman(coords="BL", M=M, a=0.5, Q=0.)
    mkn0 = KerrNewman(coords="BL", M=M, a=0., Q=0.)

    ms_mat = ms.metric_covariant(x_vec)
    mk_mat = mk.metric_covariant(x_vec)
    mk0_mat = mk0.metric_covariant(x_vec)
    mkn_mat = mkn.metric_covariant(x_vec)
    mkn0_mat = mkn0.metric_covariant(x_vec)

    v4vec_s = four_velocity(ms_mat, *v_vec, time_like)
    v4vec_k = four_velocity(mk_mat, *v_vec, time_like)
    v4vec_k0 = four_velocity(mk0_mat, *v_vec, time_like)
    v4vec_kn = four_velocity(mkn_mat, *v_vec, time_like)
    v4vec_kn0 = four_velocity(mkn0_mat, *v_vec, time_like)

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
    M = 1e24
    x_vec = np.array([1.0, np.pi / 2, 0.1])
    v_vec = np.array([-0.1, -0.01, 0.05])

    ms = Schwarzschild(M=M)
    ms_mat = ms.metric_covariant(x_vec)

    vt_s = v0(ms_mat, *v_vec)

    # Calculated by brute force
    A = ms_mat[0, 0]
    C = ms_mat[1, 1] * v_vec[0]**2 + ms_mat[2, 2] * v_vec[1]**2 + ms_mat[3, 3] * v_vec[2]**2 - _c ** 2
    D = - 4 * A * C
    vt_sb = np.sqrt(D) / (2 * A)

    assert_allclose(vt_s, vt_sb, rtol=1e-8)


def test_compare_vt_schwarzschild_kerr_kerrnewman():
    """
    Tests, whether the timelike component of 4-Velocity in KerrNewman Metric is the same as that \
    in Kerr Metric, in the limit Q -> 0 and if it becomes the same as that in Schwarzschild \
    Metric, in the limits, a -> 0 & Q -> 0

    """
    M = 1e24
    x_vec = np.array([1.0, np.pi / 2, 0.1])
    v_vec = np.array([-0.1, -0.01, 0.05])

    ms = Schwarzschild(M=M)
    mk = Kerr(coords="BL", M=M, a=0.5)
    mk0 = Kerr(coords="BL", M=M, a=0.)
    mkn = KerrNewman(coords="BL", M=M, a=0.5, Q=0.)
    mkn0 = KerrNewman(coords="BL", M=M, a=0., Q=0.)

    ms_mat = ms.metric_covariant(x_vec)
    mk_mat = mk.metric_covariant(x_vec)
    mk0_mat = mk0.metric_covariant(x_vec)
    mkn_mat = mkn.metric_covariant(x_vec)
    mkn0_mat = mkn0.metric_covariant(x_vec)

    vt_s = v0(ms_mat, *v_vec)
    vt_k = v0(mk_mat, *v_vec)
    vt_k0 = v0(mk0_mat, *v_vec)
    vt_kn = v0(mkn_mat, *v_vec)
    vt_kn0 = v0(mkn0_mat, *v_vec)

    assert_allclose(vt_s, vt_k0, rtol=1e-8)
    assert_allclose(vt_k, vt_kn, rtol=1e-8)
    assert_allclose(vt_kn0, vt_s, rtol=1e-8)
