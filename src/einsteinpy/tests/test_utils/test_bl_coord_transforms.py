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


@pytest.mark.parametrize(
    "pos_vec, a",
    [
        (np.array([10 / np.sqrt(2), 10 / np.sqrt(2), 0.0]), 0.7),
        (np.array([-732.0, 456, -90]), 9.0),
        (np.array([-732.0, -1456, 90]), 21),
        (np.array([200, -100, 0]), 0),
    ],
)
def test_cycle_pos(pos_vec, a):
    pos_vec2 = utils.CartesianToBL_pos(pos_vec, a)
    pos_vec3 = utils.BLToCartesian_pos(pos_vec2, a)
    assert_allclose(pos_vec, pos_vec3, rtol=0.0, atol=1e-5)


@pytest.mark.parametrize(
    "pos_vec, vel_vec, a",
    [
        (
            np.array([10 / np.sqrt(2), 10 / np.sqrt(2), 0.0]),
            np.array([100, -100, -43]),
            0.7,
        ),
        (np.array([-732.0, 456, -90]), np.array([2, 2, 2]), 9.0),
        (np.array([-732.0, -1456, 90]), np.array([3.4, -100, -43]), 21),
        (np.array([-10.0, -10, -20]), np.array([100, 100, 43]), 0),
        (np.array([200, -100, 0]), np.array([-500, 5000, 0]), 0),
    ],
)
def test_cycle_vel(pos_vec, vel_vec, a):
    pos_vec2 = utils.CartesianToBL_pos(pos_vec, a)
    vel_vec2 = utils.CartesianToBL_vel(pos_vec, vel_vec, a)
    # pos_vec3 = utils.BLToCartesian_pos(pos_vec2, a)
    vel_vec3 = utils.BLToCartesian_vel(pos_vec2, vel_vec2, a)
    assert_allclose(vel_vec, vel_vec3, rtol=0.0, atol=1e-5)
