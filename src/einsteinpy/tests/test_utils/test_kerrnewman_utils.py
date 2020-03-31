import numpy as np
import pytest
from astropy import units as u
from numpy.testing import assert_allclose

from einsteinpy import constant, utils
from einsteinpy.utils import kerr_utils, kerrnewman_utils

_c = constant.c.value
_G = constant.G.value
_cc = constant.coulombs_const.value


@pytest.fixture()
def test_input():
    c = _c
    G = _G
    Cc = _cc
    r = 0.1
    theta = 4 * np.pi / 5
    M = 1e23
    a = 0.99
    return c, G, Cc, r, theta, M, a


def test_compare_kerr_kerrnewman_metric_inv(test_input):
    c, G, Cc, r, theta, M, a = test_input
    # the inverse of metric for Kerr and Kerr-Newman metric should be equal when Q=0
    scr = 2 * M * G / (c ** 2)
    a_scaled = kerr_utils.scaled_spin_factor(a, M)
    m1 = kerr_utils.metric_inv(r, theta, M, a_scaled)
    m2 = kerrnewman_utils.metric_inv(r, theta, M, a_scaled, 0.0)
    assert_allclose(m1, m2, rtol=1e-10)


def test_compare_kerr_kerrnewman_dmetric_dx(test_input):
    c, G, Cc, r, theta, M, a = test_input
    # differentiation of metric for Kerr and Kerr-Newman metric should be equal when Q=0
    scr = 2 * M * G / (c ** 2)
    a_scaled = kerr_utils.scaled_spin_factor(a, M)
    m1 = kerr_utils.dmetric_dx(r, theta, M, a_scaled)
    m2 = kerrnewman_utils.dmetric_dx(r, theta, M, a_scaled, 0.0)
    assert_allclose(m1, m2, rtol=1e-10)


def test_christoffels1(test_input):
    # compare Christoffel symbols output by optimized function and by brute force
    c, G, Cc, r, theta, M, a = test_input
    Q = 1.0
    scr = 2 * M * G / (c ** 2)
    a_scaled = kerr_utils.scaled_spin_factor(a, M)
    chl1 = kerrnewman_utils.christoffels(r, theta, M, a_scaled, Q)
    # calculate by formula
    invg = kerrnewman_utils.metric_inv(r, theta, M, a_scaled, Q)
    dmdx = kerrnewman_utils.dmetric_dx(r, theta, M, a_scaled, Q)
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


def test_compare_kerr_kerrnewman_christoffels(test_input):
    # christoffel symbols for Kerr and Kerr-Newman metric should be equal when Q=0
    c, G, Cc, r, theta, M, a = test_input
    scr = 2 * M * G / (c ** 2)
    a_scaled = kerr_utils.scaled_spin_factor(a, M)
    c1 = kerr_utils.christoffels(r, theta, M, a_scaled)
    c2 = kerrnewman_utils.christoffels(r, theta, M, a_scaled, 0.0)
    assert_allclose(c1, c2, rtol=1e-8)


def test_electric_magnetic_potential_from_em_potential_vector(test_input):
    c, G, Cc, r, theta, M, a = test_input
    Q = 15.5
    vec = kerrnewman_utils.em_potential(r, theta, 0.0, Q, M)
    cmparr = np.zeros((4,), dtype=float)
    cmparr[0] = (Q / ((c ** 2) * r)) * np.sqrt(G * constant.coulombs_const.value)
    assert_allclose(cmparr, vec, rtol=1e-8)


def test_compare_kerr_kerrnewman_time_velocity():
    # time velocity for Kerr & Kerr-Newman should be same when Q=0
    pos_vec = np.array([1.0, np.pi / 2, 0.1])
    vel_vec = np.array([-0.1, -0.01, 0.05])
    mass = 1e24 * u.kg
    a01 = 0.85
    a = kerr_utils.scaled_spin_factor(a01, mass.to(u.kg).value)
    t1 = kerr_utils.kerr_time_velocity(pos_vec, vel_vec, mass, a)
    t2 = kerrnewman_utils.kerrnewman_time_velocity(pos_vec, vel_vec, mass, a, 0.0 * u.C)
    assert_allclose(t1.value, t2.value, rtol=1e-10)


def test_maxwell_tensor_covariant_for_natural_units():
    # formula for testing from https://arxiv.org/abs/gr-qc/0409025
    M = 2e22
    r = 1.5
    theta = 3 * np.pi / 5
    Q = 10.0
    a = 0.5
    m = kerrnewman_utils.maxwell_tensor_covariant(r, theta, a, Q, M, 1.0, 1.0, 1.0)
    for t in range(16):
        i = int(t / 4) % 4
        j = t % 4
        assert_allclose(0, m[i, j] + m[j, i], rtol=0.0, atol=1e-10)
    th, r2, a2 = theta, r ** 2, a ** 2
    E_r = Q * (r2 - (a * np.cos(th)) ** 2) / ((r2 + (a * np.cos(th)) ** 2) ** 2)
    E_th = (a2) * (-Q) * np.sin(2 * th) / ((r2 + (a * np.cos(th)) ** 2) ** 2)
    B_r = 2 * a * Q * (r2 + a2) * np.cos(th) / (r * ((r2 + (a * np.cos(th)) ** 2) ** 2))
    B_th = (
        a
        * Q
        * np.sin(th)
        * (r2 - (a * np.cos(th)) ** 2)
        / (r * ((r2 + (a * np.cos(th)) ** 2) ** 2))
    )
    assert_allclose(E_r, m[0, 1], rtol=1e-8)
    assert_allclose(r * E_th, m[0, 2], rtol=1e-8)
    assert_allclose(-r * np.sin(theta) * B_th / M, m[3, 1], rtol=1e-8)
    assert_allclose((r ** 2) * np.sin(theta) * B_r / M, m[3, 2], rtol=1e-8)


def test_maxwell_tensor_contravariant_for_natural_units():
    #  Theoretical background is required to write extensive tests. Right now only skew-symmetric property is being checked.
    M = 1e22
    r = 5.5
    theta = 2 * np.pi / 5
    Q = 45.0
    a = 0.7
    m = kerrnewman_utils.maxwell_tensor_contravariant(r, theta, a, Q, M, 1.0, 1.0, 1.0)
    for t in range(16):
        i = int(t / 4) % 4
        j = t % 4
        assert_allclose(0, m[i, j] + m[j, i], rtol=0.0, atol=1e-10)
