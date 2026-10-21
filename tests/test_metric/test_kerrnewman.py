import astropy.units as u
import numpy as np
import pytest

from numpy.testing import assert_allclose
from scipy.integrate import solve_ivp

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
        t=0. * u.s,
        r=0.1 * u.m,
        theta=4 * np.pi / 5 * u.rad,
        phi=0. * u.rad,
        v_r=0. * u.m / u.s,
        v_th=0. * u.rad / u.s,
        v_p=0. * u.rad / u.s
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


def test_f_vec_matches_kerr_when_Q_and_q_zero():
    """
    Validation: KerrNewman(Q=0, q=0) must reduce to Kerr.

    f_vec should be identical for both metrics when there is no charge.
    """
    M = 6.73317655e26 * u.kg
    a = 0.2 * u.one
    bl = BoyerLindquistDifferential(
        t=0.0 * u.s,
        r=1e6 * u.m,
        theta=4 * np.pi / 5 * u.rad,
        phi=0.0 * u.rad,
        v_r=0.0 * u.m / u.s,
        v_th=0.0 * u.rad / u.s,
        v_p=2e6 * u.rad / u.s,
    )

    mk = Kerr(coords=bl, M=M, a=a)
    mkn = KerrNewman(coords=bl, M=M, a=a, Q=0.0 * u.C, q=0.0 * u.C / u.kg)

    state = np.hstack((bl.position(), bl.velocity(mk)))
    f_kerr = mk._f_vec(0.0, state)
    f_kerrnewman = mkn._f_vec(0.0, state)

    assert_allclose(f_kerr, f_kerrnewman, rtol=1e-10, err_msg="Kerr vs KerrNewman f_vec mismatch when Q=q=0")


def test_f_vec_bl_kerrnewman():
    M, a = 6.73317655e26 * u.kg, 0.2 * u.one
    Q, q = 0. * u.C, 0. * u.C / u.kg
    bl = BoyerLindquistDifferential(
        t=0. * u.s,
        r=1e6 * u.m,
        theta=4 * np.pi / 5 * u.rad,
        phi=0. * u.rad,
        v_r=0. * u.m / u.s,
        v_th=0. * u.rad / u.s,
        v_p=2e6 * u.rad / u.s
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
    r = x_vec[1]

    # Using function from module (for a = 0)
    mkn = KerrNewman(coords=bl, M=M, a=0. * u.one, Q=Q)
    mkn_pot = mkn.em_potential_covariant(x_vec)

    # Calculated using formula (for a = 0)
    calc_pot = np.zeros((4,), dtype=float)
    calc_pot[0] = (Q.value / ((_c ** 2) * r)) * np.sqrt(_G * _Cc)

    assert_allclose(mkn_pot, calc_pot, rtol=1e-8)


def test_electromagnetic_potential_contravariant(test_input):
    """
    Tests, if the calculated EM 4-Potential, in contravariant form, is the same as that \
    calculated manually

    """
    bl, M, a = test_input
    Q = 15.5 * u.C
    x_vec = bl.position()
    r, theta = x_vec[1], x_vec[2]

    mkn = KerrNewman(coords=bl, M=M, a=a, Q=Q)
    mkn_contra_mat = mkn.metric_contravariant(x_vec)

    # Using function from module
    mkn_pot_contra = mkn.em_potential_contravariant(x_vec)

    # Calculated using formula
    alpha = mkn.alpha(M, a)
    rho2 = mkn.rho(r, theta, M, a) ** 2
    r_Q = np.sqrt((Q.value ** 2 * _G * _Cc) / _c ** 4)
    fac = r * r_Q / rho2
    calc_pot_cov = np.array([fac, 0., 0., -alpha * np.sin(theta)**2 * fac], dtype=float)

    calc_pot_contra = mkn_contra_mat @ calc_pot_cov

    assert_allclose(mkn_pot_contra, calc_pot_contra, rtol=1e-8)


def test_em_tensor_covariant():
    """
    Tests, if the calculated Maxwell Tensor is the same as that from the formula
    Formula for testing from https://arxiv.org/abs/gr-qc/0409025

    """
    M, a, Q = 2e22 * u.kg, 0.5 * u.one, 10.0 * u.C

    bl = BoyerLindquistDifferential(
        t=0. * u.s,
        r=1.5 * u.m,
        theta=3 * np.pi / 5 * u.rad,
        phi=0. * u.rad,
        v_r=0. * u.m / u.s,
        v_th=0. * u.rad / u.s,
        v_p=0. * u.rad / u.s
    )

    x_vec = bl.position()
    r, theta = x_vec[1], x_vec[2]

    alpha = BaseMetric.alpha(M, a)

    # Using function from module
    mkn = KerrNewman(coords=bl, M=M, a=a, Q=Q)
    mkn_em_cov = mkn.em_tensor_covariant(x_vec)

    # Checking Skew-Symmetry of Covariant Maxwell Tensor
    assert_allclose(0., mkn_em_cov + np.transpose(mkn_em_cov), atol=1e-8)

    th, r2, alpha2 = theta, r ** 2, alpha ** 2

    # Unit Scaling factor for Charge
    scale_Q = np.sqrt(_G * _Cc) / _c ** 2

    # Electric and Magnetic Fields
    D_r = (scale_Q * Q.value * (r2 - (alpha * np.cos(th)) ** 2)) / ((r2 + (alpha * np.cos(th)) ** 2) ** 2)
    D_th = ((alpha2) * (-(scale_Q * Q.value)) * np.sin(2 * th)) / ((r2 + (alpha * np.cos(th)) ** 2) ** 2)
    H_r = (2 * alpha * (scale_Q * Q.value) * (r2 + alpha2) * np.cos(th)) / (r * ((r2 + (alpha * np.cos(th)) ** 2) ** 2))
    H_th = (
        (alpha * (scale_Q * Q.value) * np.sin(th) * (r2 - (alpha * np.cos(th)) ** 2)) /
        (r * ((r2 + (alpha * np.cos(th)) ** 2) ** 2))
    )

    assert_allclose(D_r, mkn_em_cov[0, 1], rtol=1e-8)
    assert_allclose(r * D_th, mkn_em_cov[0, 2], rtol=1e-8)
    assert_allclose(r * np.sin(theta) * H_th, mkn_em_cov[1, 3], rtol=1e-8)
    assert_allclose(-(r ** 2) * np.sin(theta) * H_r, mkn_em_cov[2, 3], rtol=1e-8)


def test_em_tensor_contravariant():
    """
    Tests skew-symmetric property of EM Tensor

    Theoretical background is required to write extensive tests. \
    Right now only skew-symmetric property is being checked.

    """
    M, a, Q = 1e22 * u.kg, 0.7 * u.one, 45.0 * u.C

    bl = BoyerLindquistDifferential(
        t=0. * u.s,
        r=5.5 * u.m,
        theta=2 * np.pi / 5 * u.rad,
        phi=0. * u.rad,
        v_r=200. * u.m / u.s,
        v_th=9. * u.rad / u.s,
        v_p=10. * u.rad / u.s
    )

    x_vec = bl.position()

    # Using function from module
    mkn = KerrNewman(coords=bl, M=M, a=a, Q=Q)
    mkn_em_contra = mkn.em_tensor_contravariant(x_vec)

    assert_allclose(0., mkn_em_contra + np.transpose(mkn_em_contra), atol=1e-8)


def _kerrnewman_geometric(a, Q_M, q_m):
    """
    KerrNewman for a solar mass hole, with the hole charge Q_M given in units
    of M and the test particle charge to mass ratio q_m given in geometrized
    units. Returns the metric and the gravitational radius GM / c^2.
    """
    M = 1.989e30
    rg = _G * M / _c**2
    Q = Q_M * rg * _c**2 / np.sqrt(_G * _Cc)
    q = q_m * np.sqrt(_G / _Cc)
    bl = BoyerLindquistDifferential(
        t=0.0 * u.s,
        r=10 * rg * u.m,
        theta=np.pi / 2 * u.rad,
        phi=0.0 * u.rad,
        v_r=0.0 * u.m / u.s,
        v_th=0.0 * u.rad / u.s,
        v_p=0.0 * u.rad / u.s,
    )
    mkn = KerrNewman(coords=bl, M=M * u.kg, a=a * u.one, Q=Q * u.C, q=q * u.C / u.kg)
    return mkn, rg


def test_metric_includes_charge():
    """
    g_tt of Reissner-Nordstrom is c^2 (1 - 2M / r + Q^2 / r^2)
    """
    mkn, rg = _kerrnewman_geometric(a=0.0, Q_M=0.6, q_m=0.0)
    x_vec = np.array([0.0, 3 * rg, np.pi / 2, 0.0])

    g_tt = mkn.metric_covariant(x_vec)[0, 0]

    assert_allclose(g_tt / _c**2, 1 - 2 / 3 + 0.36 / 9, rtol=1e-12)


def test_em_tensor_is_exterior_derivative_of_potential():
    """
    F_mu_nu = d_mu A_nu - d_nu A_mu, checked with central differences
    """
    mkn, rg = _kerrnewman_geometric(a=0.6, Q_M=0.5, q_m=0.0)
    x_vec = np.array([0.0, 4 * rg, 1.0, 0.0])

    dA = np.zeros((4, 4))
    for mu, h in ((1, 1e-6 * x_vec[1]), (2, 1e-6)):
        xp, xm = x_vec.copy(), x_vec.copy()
        xp[mu] += h
        xm[mu] -= h
        dA[mu] = (mkn.em_potential_covariant(xp) - mkn.em_potential_covariant(xm)) / (
            2 * h
        )

    F = mkn.em_tensor_covariant(x_vec)

    assert_allclose(F, dA - dA.T, rtol=1e-6, atol=1e-12 * np.abs(F).max())


@pytest.mark.parametrize("Q_M, r_M", [(1.0, 1.5), (1.0, 10.0), (0.6, 3.0)])
def test_static_charged_particle_equilibrium_reissner_nordstrom(Q_M, r_M):
    """
    A charged particle at rest in Reissner-Nordstrom is in equilibrium when
    q / m = (M - Q^2 / r) / (Q sqrt(f)), with f = 1 - 2M / r + Q^2 / r^2,
    in geometrized units. For Q = M, this is q / m = 1 at every r
    (Majumdar-Papapetrou).

    """
    f = 1 - 2 / r_M + Q_M**2 / r_M**2
    q_m = (1 - Q_M**2 / r_M) / (Q_M * np.sqrt(f))
    mkn, rg = _kerrnewman_geometric(a=0.0, Q_M=Q_M, q_m=q_m)
    x_vec = np.array([0.0, r_M * rg, np.pi / 2, 0.0])
    u_t = _c / np.sqrt(mkn.metric_covariant(x_vec)[0, 0])
    state = np.hstack((x_vec, [u_t, 0.0, 0.0, 0.0]))

    f_vec = mkn._f_vec(0.0, state)
    gravity = mkn.christoffels(x_vec)[1, 0, 0] * u_t**2

    assert_allclose(f_vec[5] / gravity, 0.0, atol=1e-10)


def test_conserved_quantities_charged_orbit():
    """
    E = pi_t and L = pi_phi, with pi_mu = u_mu + (q / m) A_mu, are conserved
    along a bound charged orbit in strong field Kerr-Newman.

    """
    mkn, rg = _kerrnewman_geometric(a=0.9, Q_M=0.4, q_m=0.5)
    r0 = 8 * rg
    x0 = np.array([0.0, r0, np.pi / 2 - 0.3, 0.0])
    g0 = mkn.metric_covariant(x0)
    v_r, v_th, v_p = 0.03 * _c, 0.03 * _c / r0, 0.36 * _c / r0
    A = g0[0, 0]
    B = 2 * g0[0, 3] * v_p
    C = g0[1, 1] * v_r**2 + g0[2, 2] * v_th**2 + g0[3, 3] * v_p**2 - _c**2
    u_t = (-B + np.sqrt(B**2 - 4 * A * C)) / (2 * A)
    state0 = np.hstack((x0, [u_t, v_r, v_th, v_p]))

    # em_potential_covariant is in geometrized units with ct as x^0
    to_SI = _c * np.sqrt(_Cc / _G) * np.array([_c, 1.0, 1.0, 1.0])

    def pi_t_phi(state):
        x_vec, u_vec = state[:4], state[4:]
        pi = mkn.metric_covariant(x_vec) @ u_vec
        pi += mkn.q.value * to_SI * mkn.em_potential_covariant(x_vec)
        return pi[[0, 3]]

    T = 1000 * rg / _c
    sol = solve_ivp(mkn._f_vec, (0, T), state0, method="DOP853", rtol=1e-10, atol=1e-6)
    assert sol.success

    pi0 = pi_t_phi(state0)
    for state in sol.y.T:
        assert_allclose(pi_t_phi(state), pi0, rtol=1e-8)
