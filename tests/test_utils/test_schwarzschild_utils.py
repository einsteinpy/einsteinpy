import numpy as np
import pytest
from astropy import units as u
from numpy.testing import assert_allclose

from einsteinpy import constant, utils
from einsteinpy.utils import kerr_utils, schwarzschild_utils

_c = constant.c.value
_G = constant.G.value


def test_compare_kerr_and_schwarzschild_metric():
    r = 99.9
    theta = 5 * np.pi / 6
    M = 6.73317655e26
    # Kerr metric would reduce to Schwarzschild metric under limits a=0
    mk = kerr_utils.metric(r, theta, M, 0.0)
    ms = schwarzschild_utils.metric(r, theta, M)
    assert_allclose(mk, ms, rtol=1e-8)
