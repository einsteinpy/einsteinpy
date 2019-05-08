import numpy as np
import pytest
from astropy import units as u
from einsteinpy import constant, utils
from einsteinpy.utils import kerr_utils, schwarzschild_utils
from numpy.testing import assert_allclose

_c = constant.c.value
_G = constant.G.value


@pytest.mark.parametrize("r, theta, Rs", [(99.9, 5 * np.pi / 6, 9e-3)])
def test_compare_kerr_and_schwarzschild_metric(r, theta, Rs):
    # Kerr metric would reduce to Schwarzschild metric under limits a=0
    mk = kerr_utils.metric(_c, r, theta, Rs, 0.0)
    ms = schwarzschild_utils.metric(_c, r, theta, Rs)
    assert_allclose(mk, ms, rtol=1e-8)
