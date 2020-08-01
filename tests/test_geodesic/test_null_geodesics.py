import warnings

import numpy as np
import pytest
from numpy.testing import assert_allclose

from einsteinpy.geodesic import Nulllike
from einsteinpy.geodesic.utils import _g_dd


def test_str_repr_match():
    pos = np.array([0, 30., np.pi / 2, np.pi / 2])
    vel = np.array([-0.1, 0., 0.002])
    a = 0.9
    end_lambda = 200
    max_steps = 200
    return_cartesian = True

    geod = Nulllike(
        pos,
        vel,
        a,
        end_lambda,
        max_steps,
        return_cartesian
    )

    assert str(geod) == repr(geod)


@pytest.fixture()
def dummy_data():
    pos = np.array([0, 30., np.pi / 2, np.pi / 2])
    vel = np.array([-0.2, 0., 0.002])
    a = 0.9
    end_lambda = 200
    max_steps = 200
    return_cartesian = True

    return pos, vel, a, end_lambda, max_steps, return_cartesian


def test_NullGeodesic_has_state(dummy_data):
    pos, vel, a, end_lambda, max_steps, return_cartesian = dummy_data
    geod = Nulllike(
        pos,
        vel,
        a,
        end_lambda,
        max_steps,
        return_cartesian
    )

    assert isinstance(geod.state, np.ndarray)


def test_NullGeodesic_has_trajectory(dummy_data):
    pos, vel, a, end_lambda, max_steps, return_cartesian = dummy_data
    geod = Nulllike(
        pos,
        vel,
        a,
        end_lambda,
        max_steps,
        return_cartesian
    )

    assert isinstance(geod.trajectory[1], np.ndarray)


def test_norm_u_is_0(dummy_data):
    pos = np.array([0, 30., np.pi / 2, np.pi / 2])
    vel = np.array([-0.2, 0., 0.002])
    a = 0.9
    end_lambda = 200
    max_steps = 200
    return_cartesian = False

    geod = Nulllike(
        pos,
        vel,
        a,
        end_lambda,
        max_steps,
        return_cartesian
    )

    arr = geod.trajectory[1]
    a = geod.a
    for y in arr:
        u = y[4:]
        g_dd = _g_dd(y[1], y[2], a)

    # 1e-2 is too high, but we have to go with it for now
    assert_allclose(u @ g_dd @ u, 0., atol=1e-2)


def test_invalid_a_raises_ValueError():
    pos = np.array([0, 30., np.pi / 2, np.pi / 2])
    vel = np.array([-0.2, 0., 0.002])
    a = 1.5

    try:
        geod = Nulllike(
            position=pos,
            velocity=vel,
            a=a
        )
        assert False
    except ValueError:
        assert True
