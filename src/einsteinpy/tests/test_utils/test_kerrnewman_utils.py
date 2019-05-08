import numpy as np
import pytest
from astropy import units as u
from numpy.testing import assert_allclose

from einsteinpy import constant, utils
from einsteinpy.utils import kerr_utils, kerrnewman_utils

_c = constant.c.value
_G = constant.G.value
_cc = constant.coulombs_const.value

cos, sin = np.cos, np.sin


@pytest.mark.parametrize(
    "c, G, Cc, r, theta, M, a", [(_c, _G, _cc, 0.1, 4 * np.pi / 5, 1e23, 0.99)]
)
def test_compare_kerr_kerrnewman_metric_inv(c, G, Cc, r, theta, M, a):
    # inverse of metric for kerr and kerr-newman metric should be equal when Q=0
    scr = M * G / (c ** 2)
    a_scaled = kerr_utils.scaled_spin_factor(a, M, c, G)
    m1 = kerr_utils.metric_inv(c, r, theta, scr, a_scaled)
    m2 = kerrnewman_utils.metric_inv(c, G, Cc, r, theta, scr, a_scaled, 0.0)
    assert_allclose(m1, m2, rtol=1e-10)


@pytest.mark.parametrize(
    "c, G, Cc, r, theta, M, a", [(_c, _G, _cc, 0.1, 4 * np.pi / 5, 1e23, 0.99)]
)
def test_compare_kerr_kerrnewman_dmetric_dx(c, G, Cc, r, theta, M, a):
    # differentiation of metric for kerr and kerr-newman metric should be equal when Q=0
    scr = M * G / (c ** 2)
    a_scaled = kerr_utils.scaled_spin_factor(a, M, c, G)
    m1 = kerr_utils.dmetric_dx(c, r, theta, scr, a_scaled)
    m2 = kerrnewman_utils.dmetric_dx(c, G, Cc, r, theta, scr, a_scaled, 0.0)
    assert_allclose(m1, m2, rtol=1e-10)


@pytest.mark.parametrize(
    "c, G, Cc, r, theta, M, a, Q", [(_c, _G, _cc, 0.1, 4 * np.pi / 5, 1e23, 0.99, 1.0)]
)
def test_christoffels1(c, G, Cc, r, theta, M, a, Q):
    # compare christoffel symbols output by optimized function and by brute force
    scr = M * G / (c ** 2)
    a_scaled = kerr_utils.scaled_spin_factor(a, M, c, G)
    chl1 = kerrnewman_utils.christoffels(c, G, Cc, r, theta, scr, a_scaled, Q)
    # calculate by formula
    invg = kerrnewman_utils.metric_inv(c, G, Cc, r, theta, scr, a_scaled, Q)
    dmdx = kerrnewman_utils.dmetric_dx(c, G, Cc, r, theta, scr, a_scaled, Q)
    chl2 = np.zeros(shape=(4, 4, 4), dtype=float)
    tmp = np.array([i for i in range(4 ** 3)])
    for t in tmp:
        i = int(t / (4 ** 2)) % 4
        k = int(t / 4) % 4
        l = t % 4
        for m in range(4):
            chl2[i, k, l] += invg[i, m] * (
                dmdx[l, m, k] + dmdx[k, m, l] - dmdx[m, k, l]
            )
    chl2 = np.multiply(chl2, 0.5)
    assert_allclose(chl2, chl1, rtol=1e-10)


@pytest.mark.parametrize(
    "c, G, Cc, r, theta, M, a", [(_c, _G, _cc, 0.1, 4 * np.pi / 5, 1e23, 0.99)]
)
def test_compare_kerr_kerrnewman_christoffels(c, G, Cc, r, theta, M, a):
    # christoffel symbols for kerr and kerr-newman metric should be equal when Q=0
    scr = M * G / (c ** 2)
    a_scaled = kerr_utils.scaled_spin_factor(a, M, c, G)
    c1 = kerr_utils.christoffels(c, r, theta, scr, a_scaled)
    c2 = kerrnewman_utils.christoffels(c, G, Cc, r, theta, scr, a_scaled, 0.0)
    assert_allclose(c1, c2, rtol=1e-8)


@pytest.mark.parametrize(
    "c, G, Cc, r, theta, M, Q", [(_c, _G, _cc, 0.1, 4 * np.pi / 5, 1e23, 15.5)]
)
def test_electric_magnetic_potential_from_em_potential_vector(c, G, Cc, r, theta, M, Q):
    vec = kerrnewman_utils.em_potential(c, G, Cc, r, theta, 0.0, Q, M)
    cmparr = np.zeros((4,), dtype=float)
    cmparr[0] = (Q / ((c ** 2) * r)) * np.sqrt(G * constant.coulombs_const.value)
    assert_allclose(cmparr, vec, rtol=1e-8)


@pytest.mark.parametrize(
    "pos_vec, vel_vec, mass, a01",
    [
        (
            np.array([1.0, np.pi / 2, 0.1]),
            np.array([-0.1, -0.01, 0.05]),
            1e24 * u.kg,
            0.85,
        )
    ],
)
def test_compare_kerr_kerrnewman_time_velocity(pos_vec, vel_vec, mass, a01):
    # time velocity for kerr & kerr-newman should be same when Q=0
    a = kerr_utils.scaled_spin_factor(
        a01, mass.to(u.kg).value, constant.c.value, constant.G.value
    )
    t1 = kerr_utils.kerr_time_velocity(pos_vec, vel_vec, mass, a)
    t2 = kerrnewman_utils.kerrnewman_time_velocity(pos_vec, vel_vec, mass, a, 0.0 * u.C)
    assert_allclose(t1.value, t2.value, rtol=1e-10)


@pytest.mark.parametrize("M, r, theta, Q, a", [(2e22, 1.5, 3 * np.pi / 5, 10.0, 0.5)])
def test_maxwell_tensor_covariant_for_natural_units(M, r, theta, Q, a):
    # formula for testing from https://arxiv.org/abs/gr-qc/0409025
    m = kerrnewman_utils.maxwell_tensor_covariant(1.0, 1.0, 1.0, r, theta, a, Q, M)
    for t in range(16):
        i = int(t / 4) % 4
        j = t % 4
        assert_allclose(0, m[i, j] + m[j, i], rtol=0.0, atol=1e-10)
    th, r2, a2 = theta, r ** 2, a ** 2
    E_r = Q * (r2 - (a * cos(th)) ** 2) / ((r2 + (a * cos(th)) ** 2) ** 2)
    E_th = (a2) * (-Q) * sin(2 * th) / ((r2 + (a * cos(th)) ** 2) ** 2)
    B_r = 2 * a * Q * (r2 + a2) * cos(th) / (r * ((r2 + (a * cos(th)) ** 2) ** 2))
    B_th = (
        a
        * Q
        * sin(th)
        * (r2 - (a * cos(th)) ** 2)
        / (r * ((r2 + (a * cos(th)) ** 2) ** 2))
    )
    assert_allclose(E_r, m[0, 1], rtol=1e-8)
    assert_allclose(r * E_th, m[0, 2], rtol=1e-8)
    assert_allclose(-r * np.sin(theta) * B_th / M, m[3, 1], rtol=1e-8)
    assert_allclose((r ** 2) * np.sin(theta) * B_r / M, m[3, 2], rtol=1e-8)


@pytest.mark.parametrize("M, r, theta, Q, a", [(1e22, 5.5, 2 * np.pi / 5, 45.0, 0.7)])
def test_maxwell_tensor_contravariant_for_natural_units(M, r, theta, Q, a):
    # Theoritical background required to write extensive test. Right now only skew-symettric property is being checked.
    m = kerrnewman_utils.maxwell_tensor_contravariant(1.0, 1.0, 1.0, r, theta, a, Q, M)
    for t in range(16):
        i = int(t / 4) % 4
        j = t % 4
        assert_allclose(0, m[i, j] + m[j, i], rtol=0.0, atol=1e-10)
