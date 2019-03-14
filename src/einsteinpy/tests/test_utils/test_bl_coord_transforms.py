import numpy as np
import pytest
from astropy import units as u
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


@pytest.mark.parametrize(
    "pos_vec, vel_vec, a, ans_vel_vec",
    [
        (
            [20 * u.m, 90 * u.deg, 45 * u.deg],
            [0.0 * u.km / u.s, 0.0 * u.rad / u.s, 10 * u.rad / u.s],
            0.0,
            np.array([-200 / np.sqrt(2), 200 / np.sqrt(2), 0]),
        )
    ],
)
def test_BL2C_units(pos_vec, vel_vec, a, ans_vel_vec):
    a = utils.BL2C_units(pos_vec, vel_vec, a)
    tmp = np.array([t.value for t in a[1]])
    assert_allclose(ans_vel_vec, tmp, rtol=0.0, atol=1e-5)


@pytest.mark.parametrize(
    "vec, a",
    [
        (
            np.array(
                [
                    [21.0, 300, 0.33, 4.0, 1.0, -10.0, -1, -2],
                    [1.0, 3, 0.1, 5.0, 1.0, 100.0, -11, 2],
                    [1.0, 3, 0.1, 0, 1.0, 0, -11, 2],
                ]
            ),
            0.4,
        )
    ],
)
def test_BL2C_8dim(vec, a):
    list1 = list()
    for v in vec:
        nv = np.hstack(
            (
                v[0],
                utils.BLToCartesian_pos(v[1:4], a),
                v[4],
                utils.BLToCartesian_vel(v[1:4], v[5:8], a),
            )
        )
        list1.append(nv)
    arr1 = np.array(list1)
    arr2 = utils.BL2C_8dim(vec, a)
    assert_allclose(arr1, arr2, rtol=0.0, atol=1e-5)
