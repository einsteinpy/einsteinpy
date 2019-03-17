import collections

import numpy as np
import pytest
from astropy import units as u
from numpy.testing import assert_allclose

from einsteinpy import constant, utils
from einsteinpy.utils import kerr_utils


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
    _scr = utils.schwarzschild_radius(M * u.kg).value
    a1 = kerr_utils.event_horizon(_scr, a, np.pi / 4)
    a2 = kerr_utils.event_horizon(_scr, a, np.pi / 4, "Spherical")
    assert_allclose(a1, a2, rtol=0.0, atol=1e-5)
    assert_allclose(a1[0], _scr, rtol=0.0, atol=1e-5)


def test_ergosphere_for_nonrotating_case():
    M = 5e27
    a = 0.0
    _scr = utils.schwarzschild_radius(M * u.kg).value
    a1 = kerr_utils.ergosphere(_scr, a, np.pi / 5)
    a2 = kerr_utils.ergosphere(_scr, a, np.pi / 5, "Spherical")
    assert_allclose(a1, a2, rtol=0.0, atol=1e-5)
    assert_allclose(a1[0], _scr, rtol=0.0, atol=1e-5)


@pytest.mark.parametrize("c, r, theta, Rs, a", [(3e8, 100.0, np.pi / 5, 1.0, 0.2)])
def test_christoffels(c, r, theta, Rs, a):
    # output produced by the optimized function
    chl1 = kerr_utils.christoffels(c, r, theta, Rs, a)
    # calculate by formula
    invg = kerr_utils.metric_inv(c, r, theta, Rs, a)
    dmdx = kerr_utils.dmetric_dx(c, r, theta, Rs, a)
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
    # compare
    assert_allclose(chl2, chl1, rtol=1e-8)
