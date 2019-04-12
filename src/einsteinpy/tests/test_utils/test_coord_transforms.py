import numpy as np
import pytest
from numpy.testing import assert_allclose

from einsteinpy import utils


def test_CartesianToSpherical_pos():
    pos_vec = np.array([20.0, 311.0, 210.0])
    ans_vec = np.array([375.79382645275, 0.9778376650369, 1.5065760775947])
    ans = utils.CartesianToSpherical_pos(pos_vec)
    assert_allclose(ans, ans_vec, rtol=0.0, atol=1e-5)


def test_SphericalToCartesian_pos():
    pos_vec = np.array([12, 20 * np.pi / 180, 60 * np.pi / 180])
    ans_vec = np.array([2.052121, 3.55437, 11.27631])
    ans = utils.SphericalToCartesian_pos(pos_vec)
    assert_allclose(ans, ans_vec, rtol=0.0, atol=1e-5)


def test_SphericalToCartesian_vel():
    pos_vec = np.array([20, np.pi / 2, np.pi / 4])
    vel_vec = np.array([0.0, 0.0, 10])
    ans_vec = np.array([-200 / np.sqrt(2), 200 / np.sqrt(2), 0])
    ans = utils.SphericalToCartesian_vel(pos_vec, vel_vec)
    assert_allclose(ans, ans_vec, rtol=0.0, atol=1e-5)


def test_CartesianToSpherical_vel():
    pos_vec = np.array([10 / np.sqrt(2), 10 / np.sqrt(2), 0.0])
    vel_vec = np.array([-190 / np.sqrt(2), 210 / np.sqrt(2), 200.0])
    ans_vec = np.array([10.0, -20.0, 20.0])
    ans = utils.CartesianToSpherical_vel(pos_vec, vel_vec)
    assert_allclose(ans, ans_vec, rtol=0.0, atol=1e-5)


def test_S2C_8dim():
    vec = np.array(
        [
            [21.0, 300, 0.33, 4.0, 1.0, -10.0, -1, -2],
            [1.0, 3, 0.1, 5.0, 1.0, 100.0, -11, 2],
            [1.0, 3, 0.1, 0, 1.0, 0, -11, 2],
        ]
    )
    list1 = list()
    for v in vec:
        nv = np.hstack(
            (
                v[0],
                utils.SphericalToCartesian_pos(v[1:4]),
                v[4],
                utils.SphericalToCartesian_vel(v[1:4], v[5:8]),
            )
        )
        list1.append(nv)
    arr1 = np.array(list1)
    arr2 = utils.S2C_8dim(vec)
    assert_allclose(arr1, arr2, rtol=0.0, atol=1e-5)
