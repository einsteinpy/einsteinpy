import collections

import numpy as np
import pytest
from astropy import units as u
from numpy.testing import assert_allclose

from einsteinpy import constant
from einsteinpy.utils import kerr_utils, schwarzschild_radius_dimensionless


def test_nonzero_christoffels():
    l1 = kerr_utils.nonzero_christoffels()
    l2 = kerr_utils.nonzero_christoffels_list
    assert collections.Counter(l1) == collections.Counter(l2)


def test_spin_factor():
    J, M, c = 9e10, 5e24, constant.c.value
    assert_allclose(J / (M * c), kerr_utils.spin_factor(J, M, c), rtol=0.0, atol=1e-10)


def test_event_horizon_for_nonrotating_case():
    M = 5e27
    a = 0.0
    _scr = schwarzschild_radius_dimensionless(M)
    a1 = kerr_utils.event_horizon(M, a, np.pi / 4)
    a2 = kerr_utils.event_horizon(M, a, np.pi / 4, "Spherical")
    assert_allclose(a1, a2, rtol=0.0, atol=1e-5)
    assert_allclose(a1[0], _scr, rtol=0.0, atol=1e-5)


def test_radius_ergosphere_for_nonrotating_case():
    M = 5e27
    a = 0.0
    _scr = schwarzschild_radius_dimensionless(M)
    a1 = kerr_utils.radius_ergosphere(M, a, np.pi / 5)
    a2 = kerr_utils.radius_ergosphere(M, a, np.pi / 5, "Spherical")
    assert_allclose(a1, a2, rtol=0.0, atol=1e-5)
    assert_allclose(a1[0], _scr, rtol=0.0, atol=1e-5)


def test_christoffels():
    # output produced by the optimized function
    c = 3e8
    r = 100.0
    theta = np.pi / 5
    M = 6.73317655e26
    a = 0.2
    chl1 = kerr_utils.christoffels(r, theta, M, a, c)
    # calculate by formula
    invg = kerr_utils.metric_inv(r, theta, M, a, c)
    dmdx = kerr_utils.dmetric_dx(r, theta, M, a, c)
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
    assert_allclose(chl2, chl1, rtol=1e-8)


@pytest.mark.parametrize("a", [1.1, -0.1])
def test_scaled_spin_factor_raises_error(a):
    try:
        kerr_utils.scaled_spin_factor(a, 3e20)
        assert False
    except ValueError:
        assert True


def test_scaled_spin_factor():
    a = 0.8
    M = 5e24
    a1 = kerr_utils.scaled_spin_factor(a, M)
    a2 = schwarzschild_radius_dimensionless(M) * a * 0.5
    assert_allclose(a2, a1, rtol=1e-9)
