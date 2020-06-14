import numpy as np
import pytest

from numpy.testing import assert_allclose

from einsteinpy import constant
from einsteinpy.metric import GenericMetric, Kerr, KerrNewman

_c = constant.c.value
_G = constant.G.value
_Cc = constant.coulombs_const.value


@pytest.fixture()
def test_input():
    r = 0.1
    theta = 4 * np.pi / 5
    M = 1e23
    a = 0.99

    return r, theta, M, a


def test_compare_kerr_kerrnewman_dmetric_dx(test_input):
    """
    Tests, if the metric derivatives for Kerr & Kerr-Newman metrics match, when Q -> 0

    """
    r, theta, M, a = test_input
    x_vec = np.array([0., r, theta, 0.])

    mk = Kerr(coords="BL", M=M, a=a)
    mkn = KerrNewman(coords="BL", M=M, a=a, Q=0.)
    mkdx = mk._dg_dx_bl(x_vec)
    mkndx = mkn._dg_dx_bl(x_vec)

    assert_allclose(mkdx, mkndx, rtol=1e-10)


def test_christoffels_kerr_newman(test_input):
    """
    Compares output produced by optimized function, with that, produced via general method (formula)

    """
    r, theta, M, a = test_input
    Q = 1.0
    x_vec = np.array([0., r, theta, 0.])

    # Output produced by the optimized function
    mkn = KerrNewman(coords="BL", M=M, a=a, Q=Q)
    chl1 = mkn.christoffels(x_vec)

    # Calculated using formula
    g_contra = mkn.metric_contravariant(x_vec)
    dgdx = mkn._dg_dx_bl(x_vec)
    chl2 = np.zeros(shape=(4, 4, 4), dtype=float)
    tmp = np.array([i for i in range(4 ** 3)])
    for t in tmp:
        i = int(t / (4 ** 2)) % 4
        k = int(t / 4) % 4
        l = t % 4
        for m in range(4):
            chl2[i, k, l] += g_contra[i, m] * (
                dgdx[l, m, k] + dgdx[k, m, l] - dgdx[m, k, l]
            )
    chl2 = np.multiply(chl2, 0.5)
    
    assert_allclose(chl2, chl1, rtol=1e-10)


def test_compare_kerr_kerrnewman_christoffels(test_input):
    """
    Compares KerrNewman Christoffel Symbols, with that of Kerr metric, when Q -> 0

    """
    r, theta, M, a = test_input
    x_vec = np.array([0., r, theta, 0.])

    mk = Kerr(coords="BL", M=M, a=a)
    mkn = KerrNewman(coords="BL", M=M, a=a, Q=0.)
    mk_chl = mk.christoffels(x_vec)
    mkn_chl = mkn.christoffels(x_vec)

    assert_allclose(mk_chl, mkn_chl, rtol=1e-8)


def test_electromagnetic_potential_from_em_potential_vector(test_input):
    """
    Tests, if the calculated EM 4-Potential is the same as that from the formula

    """
    r, theta, M, a = test_input
    Q = 15.5

    # Using function from module
    mkn = KerrNewman(coords="BL", M=M, a=a, Q=Q)
    mkn_pot = mkn.em_potential_covariant(r, theta, M=M, a=0., Q=Q)

    # Calculated using formula
    calc_pot = np.zeros((4,), dtype=float)
    calc_pot[0] = (Q / ((_c ** 2) * r)) * np.sqrt(_G * _Cc)

    assert_allclose(mkn_pot, calc_pot, rtol=1e-8)


def test_em_tensor_covariant():
    """
    Tests, if the calculated Maxwell Tensor is the same as that from the formula
    Formula for testing from https://arxiv.org/abs/gr-qc/0409025

    """
    r, theta = 1.5, 3 * np.pi / 5
    M, a, Q = 2e22, 0.5, 10.0
    alpha = GenericMetric.alpha(M, a)

    # Using function from module
    mkn = KerrNewman(coords="BL", M=M, a=a, Q=Q)
    mkn_em_cov = mkn.em_tensor_covariant(r, theta, M, a, Q)

    # Checking Skew-Symmetry of Covariant Maxwell Tensor
    assert_allclose(0., mkn_em_cov + np.transpose(mkn_em_cov), atol=1e-8)
    
    th, r2, alpha2 = theta, r ** 2, alpha ** 2
    
    # Unit Scaling factor for Charge
    scale_Q = np.sqrt(_G * _Cc) / _c ** 2
    
    # Electric and Magnetic Fields
    D_r = (scale_Q * Q * (r2 - (alpha * np.cos(th)) ** 2)) / ((r2 + (alpha * np.cos(th)) ** 2) ** 2)
    D_th = ((alpha2) * (-(scale_Q * Q)) * np.sin(2 * th)) / ((r2 + (alpha * np.cos(th)) ** 2) ** 2)
    H_r = (2 * alpha * (scale_Q * Q) * (r2 + alpha2) * np.cos(th)) / (r * ((r2 + (alpha * np.cos(th)) ** 2) ** 2))
    H_th = (alpha * (scale_Q * Q) * np.sin(th) * (r2 - (alpha * np.cos(th)) ** 2)) / (r * ((r2 + (alpha * np.cos(th)) ** 2) ** 2))
    
    assert_allclose(D_r, mkn_em_cov[0, 1], rtol=1e-8)
    assert_allclose(r * D_th, mkn_em_cov[0, 2], rtol=1e-8)
    assert_allclose(-r * np.sin(theta) * H_th, mkn_em_cov[1, 3], rtol=1e-8)
    assert_allclose((r ** 2) * np.sin(theta) * H_r, mkn_em_cov[2, 3], rtol=1e-8)


def test_em_tensor_contravariant():
    """
    Tests skew-symmetric property of EM Tensor

    Theoretical background is required to write extensive tests. \
    Right now only skew-symmetric property is being checked.

    """
    #  
    r, theta = 5.5, 2 * np.pi / 5
    M, a, Q = 1e22, 0.7, 45.0

    # Using function from module
    mkn = KerrNewman(coords="BL", M=M, a=a, Q=Q)
    mkn_em_contra = mkn.em_tensor_contravariant(r, theta, M, a, Q)
    
    assert_allclose(0., mkn_em_contra + np.transpose(mkn_em_contra), atol=1e-8)
