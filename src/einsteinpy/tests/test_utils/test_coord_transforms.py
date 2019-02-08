import numpy as np
import pytest
from numpy.testing import assert_allclose

from einsteinpy import utils

@pytest.mark.parametrize("pos_vec, ans_vec",
    [
        (np.array([20.0,311.0,210.0]), np.array([375.79382645275, 0.9778376650369, 1.5065760775947])),
    ])
def test_CartesianToSpherical_pos(pos_vec, ans_vec):
    ans = utils.CartesianToSpherical_pos(pos_vec)
    assert_allclose(ans, ans_vec, rtol=0., atol=1e-5)


@pytest.mark.parametrize("pos_vec, ans_vec",
    [
        (np.array([50.,60.,-30.]), np.array([-2.350874012, -15.05812665	, -47.62064902])),
    ])
def test_SphericalToCartesian_pos(pos_vec, ans_vec):
    ans = utils.CartesianToSpherical_pos(pos_vec)
    assert_allclose(ans, ans_vec, rtol=0., atol=1e-5)