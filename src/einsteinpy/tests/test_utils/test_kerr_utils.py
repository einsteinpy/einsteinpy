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
