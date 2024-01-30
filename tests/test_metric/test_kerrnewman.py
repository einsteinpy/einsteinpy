import astropy.units as u
import numpy as np
import pytest

from numpy.testing import assert_allclose

from einsteinpy import constant
from einsteinpy.coordinates import BoyerLindquistDifferential
from einsteinpy.metric import BaseMetric, Kerr, KerrNewman

_c = constant.c.value
_G = constant.G.value
_Cc = constant.coulombs_const.value


@pytest.fixture()
def test_input():
    """
    Test input for some functions below
    """
    bl = BoyerLindquistDifferential(
        e0=0. * u.s,
        e1=0.1 * u.m,
        e2=4 * np.pi / 5 * u.rad,
        e3=0. * u.rad,
        u0=0. * u.m / u.s,
        u1=0. * u.rad / u.s,
        u2=0. * u.rad / u.s
    )

    M = 1e23 * u.kg
    a = 0.99 * u.one

    return bl, M, a


def test_compare_kerr_kerrnewman_dmetric_dx(test_input):
    """
    Tests, if the metric derivatives for Kerr & Kerr-Newman metrics match, when Q -> 0

    """
    bl, M, a = test_input
    x_vec = bl.position()

    mk = Kerr(coords=bl, M=M, a=a)
    mkn = KerrNewman(coords=bl, M=M, a=a, Q=0. * u.C)
    mkdx = mk._dg_dx_bl(x_vec)
    mkndx = mkn._dg_dx_bl(x_vec)

    assert_allclose(mkdx, mkndx, rtol=1e-8)


def test_christoffels_kerrnewman(test_input):
    """
    Compares output produced by optimized function, with that, produced via general method (formula)

    """
    bl, M, a = test_input
    Q = 1.0 * u.C
    x_vec = bl.position()

    # Output produced by the optimized function
    mkn = KerrNewman(coords=bl, M=M, a=a, Q=Q)
    chl1 = mkn.christoffels(x_vec)

    # Calculated using formula
    g_contra = mkn.metric_contravariant(x_vec)
    dgdx = mkn._dg_dx_bl(x_vec)
    chl2 = np.zeros(shape=(4, 4, 4), dtype=float)
    tmp = np.array([i for i in range(4 ** 3)])
    for t in tmp:
        i = int(t / (4 ** 2)) % 4
        k = int(t / 4) % 4
        index = t % 4
        for m in range(4):
            chl2[i, k, index] += g_contra[i, m] * (
                dgdx[index, m, k] + dgdx[k, m, index] - dgdx[m, k, index]
            )
    chl2 = np.multiply(chl2, 0.5)

    assert_allclose(chl2, chl1, rtol=1e-10)


def test_f_vec_bl_kerrnewman():
    M, a = 6.73317655e26 * u.kg, 0.2 * u.one
    Q, q = 0. * u.C, 0. * u.C / u.kg
    bl = BoyerLindquistDifferential(
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

    mk = KerrNewman(coords=bl, M=M, a=a, Q=Q, q=q)
    state = np.hstack((bl.position(), bl.velocity(mk)))

    f_vec = mk._f_vec(0., state)

    assert isinstance(f_vec, np.ndarray)
    assert_allclose(f_vec_expected, f_vec, rtol=1e-8)


def test_compare_kerr_kerrnewman_christoffels(test_input):
    """
    Compares KerrNewman Christoffel Symbols, with that of Kerr metric, when Q -> 0

    """
    bl, M, a = test_input
    x_vec = bl.position()

    mk = Kerr(coords=bl, M=M, a=a)
    mkn = KerrNewman(coords=bl, M=M, a=a, Q=0. * u.C)
    mk_chl = mk.christoffels(x_vec)
    mkn_chl = mkn.christoffels(x_vec)

    assert_allclose(mk_chl, mkn_chl, rtol=1e-8)


def test_electromagnetic_potential_from_em_potential_vector(test_input):
    """
    Tests, if the calculated EM 4-Potential is the same as that from the formula

    """
    bl, M, a = test_input
    Q = 15.5 * u.C

    x_vec = bl.position()
    e1 = x_vec[1]

    # Using function from module (for a = 0)
    mkn = KerrNewman(coords=bl, M=M, a=0. * u.one, Q=Q)
    mkn_pot = mkn.em_potential_covariant(x_vec)

    # Calculated using formula (for a = 0)
    calc_pot = np.zeros((4,), dtype=float)
    calc_pot[0] = (Q.value / ((_c ** 2) * e1)) * np.sqrt(_G * _Cc)

    assert_allclose(mkn_pot, calc_pot, rtol=1e-8)


def test_electromagnetic_potential_contravariant(test_input):
    """
    Tests, if the calculated EM 4-Potential, in contravariant form, is the same as that \
    calculated manually

    """
    bl, M, a = test_input
    Q = 15.5 * u.C
    x_vec = bl.position()
    e1, e2 = x_vec[1], x_vec[2]

    mkn = KerrNewman(coords=bl, M=M, a=a, Q=Q)
    mkn_contra_mat = mkn.metric_contravariant(x_vec)

    # Using function from module
    mkn_pot_contra = mkn.em_potential_contravariant(x_vec)

    # Calculated using formula
    alpha = mkn.alpha(M, a)
    rho2 = mkn.rho(e1, e2, M, a) ** 2
    r_Q = np.sqrt((Q.value ** 2 * _G * _Cc) / _c ** 4)
    fac = e1 * r_Q / rho2
    calc_pot_cov = np.array([fac, 0., 0., -alpha * np.sin(e2)**2 * fac], dtype=float)

    calc_pot_contra = mkn_contra_mat @ calc_pot_cov

    assert_allclose(mkn_pot_contra, calc_pot_contra, rtol=1e-8)


def test_em_tensor_covariant():
    """
    Tests, if the calculated Maxwell Tensor is the same as that from the formula
    Formula for testing from https://arxiv.org/abs/gr-qc/0409025

    """
    M, a, Q = 2e22 * u.kg, 0.5 * u.one, 10.0 * u.C

    bl = BoyerLindquistDifferential(
        e0=0. * u.s,
        e1=1.5 * u.m,
        e2=3 * np.pi / 5 * u.rad,
        e3=0. * u.rad,
        u0=0. * u.m / u.s,
        u1=0. * u.rad / u.s,
        u2=0. * u.rad / u.s
    )

    x_vec = bl.position()
    e1, e2 = x_vec[1], x_vec[2]

    alpha = BaseMetric.alpha(M, a)

    # Using function from module
    mkn = KerrNewman(coords=bl, M=M, a=a, Q=Q)
    mkn_em_cov = mkn.em_tensor_covariant(x_vec)

    # Checking Skew-Symmetry of Covariant Maxwell Tensor
    assert_allclose(0., mkn_em_cov + np.transpose(mkn_em_cov), atol=1e-8)

    th, r2, alpha2 = e2, e1 ** 2, alpha ** 2

    # Unit Scaling factor for Charge
    scale_Q = np.sqrt(_G * _Cc) / _c ** 2

    # Electric and Magnetic Fields
    D_r = (scale_Q * Q.value * (r2 - (alpha * np.cos(th)) ** 2)) / ((r2 + (alpha * np.cos(th)) ** 2) ** 2)
    D_th = ((alpha2) * (-(scale_Q * Q.value)) * np.sin(2 * th)) / ((r2 + (alpha * np.cos(th)) ** 2) ** 2)
    H_r = (2 * alpha * (scale_Q * Q.value) * (r2 + alpha2) * np.cos(th)) / (e1 * ((r2 + (alpha * np.cos(th)) ** 2) ** 2))
    H_th = (
        (alpha * (scale_Q * Q.value) * np.sin(th) * (r2 - (alpha * np.cos(th)) ** 2)) /
        (e1 * ((r2 + (alpha * np.cos(th)) ** 2) ** 2))
    )

    assert_allclose(D_r, mkn_em_cov[0, 1], rtol=1e-8)
    assert_allclose(e1 * D_th, mkn_em_cov[0, 2], rtol=1e-8)
    assert_allclose(-e1 * np.sin(e2) * H_th, mkn_em_cov[1, 3], rtol=1e-8)
    assert_allclose((e1 ** 2) * np.sin(e2) * H_r, mkn_em_cov[2, 3], rtol=1e-8)


def test_em_tensor_contravariant():
    """
    Tests skew-symmetric property of EM Tensor

    Theoretical background is required to write extensive tests. \
    Right now only skew-symmetric property is being checked.

    """
    M, a, Q = 1e22 * u.kg, 0.7 * u.one, 45.0 * u.C

    bl = BoyerLindquistDifferential(
        e0=0. * u.s,
        e1=5.5 * u.m,
        e2=2 * np.pi / 5 * u.rad,
        e3=0. * u.rad,
        u0=200. * u.m / u.s,
        u1=9. * u.rad / u.s,
        u2=10. * u.rad / u.s
    )

    x_vec = bl.position()

    # Using function from module
    mkn = KerrNewman(coords=bl, M=M, a=a, Q=Q)
    mkn_em_contra = mkn.em_tensor_contravariant(x_vec)

    assert_allclose(0., mkn_em_contra + np.transpose(mkn_em_contra), atol=1e-8)
