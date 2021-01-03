import warnings

import numpy as np
import pytest
from numpy.testing import assert_allclose

from einsteinpy.geodesic import Geodesic, Nulllike, Timelike


@pytest.fixture()
def dummy_timegeod():
    """
    Equatorial Capture

    """
    return Timelike(
        metric="Kerr",
        metric_params=(0.9,),
        position=[2.15, np.pi / 2, 0.],
        momentum=[0., 0., 1.5],
        steps=50,
        delta=0.5,
        omega=0.01,  # Close orbit
        return_cartesian=True,
        suppress_warnings=True,
    )


@pytest.fixture()
def dummy_nullgeod():
    """
    Equatorial Geodesic

    """
    return Nulllike(
        metric="Kerr",
        metric_params=(0.5,),
        position=[4., np.pi / 2, 0.],
        momentum=[0., 0., 2.],
        steps=50,
        delta=0.5,
        return_cartesian=False,
        suppress_warnings=True,
    )


def test_str_repr(dummy_timegeod):
    geod = dummy_timegeod

    assert str(geod) == repr(geod)


def test_NotImplementedError():
    try:
        geod = Nulllike(
            metric="Ker",
            metric_params=(0.9,),
            position=[2.5, np.pi / 2, 0.],
            momentum=[0., 0., -8.5],
        )
        assert False

    except NotImplementedError:
        assert True


def test_geodesic_attributes(dummy_timegeod):
    geod = dummy_timegeod
    traj = geod.trajectory

    assert traj
    assert isinstance(traj, tuple)
    assert isinstance(traj[0], np.ndarray)
    assert isinstance(traj[1], np.ndarray)
    assert traj[0].shape[0] == traj[1].shape[0]
    assert traj[1].shape[1] == 8


def test_constant_angular_momentum(dummy_nullgeod):
    L = dummy_nullgeod.momentum[-1]

    assert_allclose(dummy_nullgeod.trajectory[1][:, -1], L, atol=1e-4, rtol=1e-4)


def test_equatorial_geodesic(dummy_nullgeod):
    theta = dummy_nullgeod.position[2]

    assert_allclose(dummy_nullgeod.trajectory[1][:, 2], theta, atol=1e-6, rtol=1e-6)


def test_constant_rad():
    geod = Timelike(
        metric="Kerr",
        metric_params=(0.99,),
        position=[4., np.pi / 3, 0.],
        momentum=[0., 0.767851, 2.],
        return_cartesian=False,
        steps=50,
        delta=1.,
    )
    r = geod.trajectory[1][:, 1]

    assert_allclose(r, 4., atol=1e-2, rtol=1e-2)


def test_kerr0_eq_sch():
    metric_params = (0.,)
    q0 = [4., np.pi / 2, 0.]
    p0 =  [0., 0., 0.]

    k = Timelike(
        metric="Kerr",
        metric_params=metric_params,
        position=q0,
        momentum=p0,
        steps=50,
        delta=0.5,
        return_cartesian=True,
        suppress_warnings=True,
    )

    s = Timelike(
        metric="Schwarzschild",
        metric_params=metric_params,
        position=q0,
        momentum=p0,
        steps=50,
        delta=0.5,
        return_cartesian=True,
        suppress_warnings=True,
    )

    assert_allclose(k.trajectory[0], s.trajectory[0], atol=1e-6, rtol=1e-6)
    assert_allclose(k.trajectory[1], s.trajectory[1], atol=1e-6, rtol=1e-6)


def test_kerr0_eq_kn00():
    metric_params = (0.5, 0.)
    q0 = [2.5, np.pi / 2, 0.]
    p0 = [0., 0., -8.5]

    k = Timelike(
        metric="Kerr",
        metric_params=metric_params,
        position=q0,
        momentum=p0,
        steps=50,
        delta=0.5,
        return_cartesian=True,
        suppress_warnings=True,
    )

    kn = Timelike(
        metric="KerrNewman",
        metric_params=metric_params,
        position=q0,
        momentum=p0,
        steps=50,
        delta=0.5,
        return_cartesian=True,
        suppress_warnings=True,
    )

    assert_allclose(k.trajectory[0], kn.trajectory[0], atol=1e-6, rtol=1e-6)
    assert_allclose(k.trajectory[1], kn.trajectory[1], atol=1e-6, rtol=1e-6)
