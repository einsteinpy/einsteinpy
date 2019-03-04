import numpy as np
import pytest
from numpy.testing import assert_allclose

from einsteinpy import utils


@pytest.mark.parametrize(
    "pos_vec, a, ans_vec",
    [
        (
            np.array([20.0, 311.0, 210.0]),
            0.0,
            np.array([375.79382645275, 0.9778376650369, 1.5065760775947]),
        )
    ],
)
def test_CartesianToBL_pos(pos_vec, a, ans_vec):
    ans = utils.CartesianToBL_pos(pos_vec, a)
    assert_allclose(ans, ans_vec, rtol=0.0, atol=1e-5)


@pytest.mark.parametrize(
    "pos_vec, a, ans_vec",
    [
        (
            np.array([12, 20 * np.pi / 180, 60 * np.pi / 180]),
            0.0,
            np.array([2.052121, 3.55437, 11.27631]),
        )
    ],
)
def test_BLtoCartesian_pos(pos_vec, a, ans_vec):
    ans = utils.BLToCartesian_pos(pos_vec, a)
    assert_allclose(ans, ans_vec, rtol=0.0, atol=1e-5)


@pytest.mark.parametrize(
    "pos_vec, vel_vec, a, ans_vec",
    [
        (
            np.array([10 / np.sqrt(2), 10 / np.sqrt(2), 0.0]),
            np.array([-190 / np.sqrt(2), 210 / np.sqrt(2), 200.0]),
            0.0,
            np.array([10.0, -20.0, 20.0]),
        )
    ],
)
def test_CartesianToBL_vel(pos_vec, vel_vec, a, ans_vec):
    ans = utils.CartesianToBL_vel(pos_vec, vel_vec, a)
    assert_allclose(ans, ans_vec, rtol=0.0, atol=1e-5)


@pytest.mark.parametrize(
    "pos_vec, vel_vec, a, ans_vec",
    [
        (
            np.array([20, np.pi / 2, np.pi / 4]),
            np.array([0.0, 0.0, 10]),
            0.0,
            np.array([-200 / np.sqrt(2), 200 / np.sqrt(2), 0]),
        )
    ],
)
def test_BLtoCartesian_vel(pos_vec, vel_vec, a, ans_vec):
    ans = utils.BLToCartesian_vel(pos_vec, vel_vec, a)
    assert_allclose(ans, ans_vec, rtol=0.0, atol=1e-5)
