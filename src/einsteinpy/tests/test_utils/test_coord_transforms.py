import numpy as np
import pytest
from numpy.testing import assert_allclose

from einsteinpy import utils


@pytest.mark.parametrize(
    "pos_vec, ans_vec",
    [
        (
            np.array([20.0, 311.0, 210.0]),
            np.array([375.79382645275, 0.9778376650369, 1.5065760775947]),
        )
    ],
)
def test_CartesianToSpherical_pos(pos_vec, ans_vec):
    ans = utils.CartesianToSpherical_pos(pos_vec)
    assert_allclose(ans, ans_vec, rtol=0.0, atol=1e-5)


@pytest.mark.parametrize(
    "pos_vec, ans_vec",
    [
        (
            np.array([12, 20 * np.pi / 180, 60 * np.pi / 180]),
            np.array([2.052121, 3.55437, 11.27631]),
        )
    ],
)
def test_SphericalToCartesian_pos(pos_vec, ans_vec):
    ans = utils.SphericalToCartesian_pos(pos_vec)
    assert_allclose(ans, ans_vec, rtol=0.0, atol=1e-5)


@pytest.mark.parametrize(
    "pos_vec, vel_vec, ans_vec",
    [
        (
            np.array([20, np.pi / 2, np.pi / 4]),
            np.array([0.0, 0.0, 10]),
            np.array([-200 / np.sqrt(2), 200 / np.sqrt(2), 0]),
        )
    ],
)
def test_SphericalToCartesian_vel(pos_vec, vel_vec, ans_vec):
    ans = utils.SphericalToCartesian_vel(pos_vec, vel_vec)
    assert_allclose(ans, ans_vec, rtol=0.0, atol=1e-5)


@pytest.mark.parametrize(
    "pos_vec, vel_vec, ans_vec",
    [
        (
            np.array([10 / np.sqrt(2), 10 / np.sqrt(2), 0.0]),
            np.array([-190 / np.sqrt(2), 210 / np.sqrt(2), 200.0]),
            np.array([10.0, -20.0, 20.0]),
        )
    ],
)
def test_CartesianToSpherical_vel(pos_vec, vel_vec, ans_vec):
    ans = utils.CartesianToSpherical_vel(pos_vec, vel_vec)
    assert_allclose(ans, ans_vec, rtol=0.0, atol=1e-5)
