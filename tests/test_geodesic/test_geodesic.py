import warnings

import numpy as np
import pytest
from numpy.testing import assert_allclose

from einsteinpy.geodesic import Geodesic, Nulllike, Timelike


@pytest.fixture()
def dummy_time_julia():
    """
    Equatorial Unstable Circular Capture

    """
    q0 = [4, np.pi / 2, 0.]
    p0 = [0., 0., 3.99999]
    a = 0.
    end_lambda = 300.
    step_size = 0.05
    julia = True

    return q0, p0, a, end_lambda, step_size, julia


@pytest.fixture()
def dummy_null_julia():
    """
    Equatorial Reverse & Escape, High Spin

    """
    q0 = [4, np.pi / 2, 0.]
    p0 = [0., 0., -2.]
    a = 0.9
    end_lambda = 150.
    step_size = 0.005
    julia = True

    return q0, p0, a, end_lambda, step_size, julia


@pytest.fixture()
def dummy_time_python():
    """
    Equatorial Spiral Capture

    """
    q0 = [2.15, np.pi / 2, 0.]
    p0 = [0., 0., -1.5]
    a = 0.
    end_lambda = 10.
    step_size = 0.005
    julia = False

    return q0, p0, a, end_lambda, step_size, julia


@pytest.fixture()
def dummy_null_python():
    """
    Equatorial Reverse & Capture

    """
    q0 = [2.5, np.pi / 2, 0.]
    p0 = [0., 0., -8.5]
    a = 0.9
    end_lambda = 10.
    step_size = 0.005
    julia = False

    return q0, p0, a, end_lambda, step_size, julia


def test_basics_python1(dummy_time_python):
    q0, p0, a, end_lambda, step_size, julia = dummy_time_python
    geod = Timelike(
        position=q0,
        momentum=p0,
        a=a,
        end_lambda=end_lambda,
        step_size=step_size,
        return_cartesian=True,
        julia=julia
    )

    assert str(geod) == repr(geod)
    assert isinstance(geod.trajectory, tuple)
    assert isinstance(geod.trajectory[0], np.ndarray)
    assert isinstance(geod.trajectory[1], np.ndarray)
    assert isinstance(geod.trajectory[0], np.ndarray)


def test_basics_python2(dummy_time_python):
    q0, p0, a, end_lambda, step_size, julia = dummy_time_python
    geod = Timelike(
        position=q0,
        momentum=p0,
        a=a,
        end_lambda=end_lambda,
        step_size=step_size,
        return_cartesian=True,
        julia=julia
    )

    assert geod.trajectory
    assert geod.trajectory[0].shape[0] == geod.trajectory[1].shape[0]
    assert geod.trajectory[1].shape[1] == 6


def test_basics_julia1(dummy_null_julia):
    q0, p0, a, end_lambda, step_size, julia = dummy_null_julia
    geod = Nulllike(
        position=q0,
        momentum=p0,
        a=a,
        end_lambda=end_lambda,
        step_size=step_size,
        return_cartesian=True,
        julia=julia
    )

    assert str(geod) == repr(geod)
    assert isinstance(geod.trajectory, tuple)
    assert isinstance(geod.trajectory[0], np.ndarray)
    assert isinstance(geod.trajectory[1], np.ndarray)
    assert isinstance(geod.trajectory[0], np.ndarray)


def test_basics_julia2(dummy_null_julia):
    q0, p0, a, end_lambda, step_size, julia = dummy_null_julia
    geod = Nulllike(
        position=q0,
        momentum=p0,
        a=a,
        end_lambda=end_lambda,
        step_size=step_size,
        return_cartesian=True,
        julia=julia
    )

    assert geod.trajectory
    assert geod.trajectory[0].shape[0] == geod.trajectory[1].shape[0]
    assert geod.trajectory[1].shape[1] == 6


def test_runtime_warning_python(dummy_time_python):
    q0, p0, a, end_lambda, step_size, julia = dummy_time_python

    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")

        geod = Geodesic(
            position=q0,
            momentum=p0,
            a=a,
            end_lambda=end_lambda,
            step_size=step_size,
            time_like=True,
            return_cartesian=True,
            julia=julia
        )

        assert len(w) == 1  # 1 warning to be shown
        assert issubclass(w[-1].category, RuntimeWarning)


def test_runtime_warning_julia(dummy_time_julia):
    q0, p0, a, end_lambda, step_size, julia = dummy_time_julia

    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")

        geod = Timelike(
            position=q0,
            momentum=p0,
            a=a,
            end_lambda=end_lambda,
            step_size=step_size,
            return_cartesian=True,
            julia=julia
        )

        assert len(w) == 1  # 1 warning to be shown
        assert issubclass(w[-1].category, RuntimeWarning)


def test_python_use_warning(dummy_null_python):
    q0, p0, a, end_lambda, step_size, julia = dummy_null_python

    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")

        geod = Nulllike(
            position=q0,
            momentum=p0,
            a=a,
            end_lambda=end_lambda,
            step_size=step_size,
            return_cartesian=True,
            julia=julia
        )

        assert len(w) == 2  # 2 warnings to be shown (as capture geodesic)
        assert issubclass(w[-2].category, RuntimeWarning)
        assert issubclass(w[-1].category, RuntimeWarning)


def test_constant_angular_momentum(dummy_null_python):
    q0, p0, a, end_lambda, step_size, julia = dummy_null_python

    geod = Nulllike(
        position=q0,
        momentum=p0,
        a=a,
        end_lambda=end_lambda,
        step_size=step_size,
        return_cartesian=True,
        julia=julia
    )

    L = p0[-1]

    assert_allclose(geod.trajectory[1][:, 5], L, atol=1e-4, rtol=1e-4)


def test_equatorial_geodesic(dummy_time_python):
    q0, p0, a, end_lambda, step_size, julia = dummy_time_python

    geod = Timelike(
        position=q0,
        momentum=p0,
        a=a,
        end_lambda=end_lambda,
        step_size=step_size,
        return_cartesian=False,
        julia=julia
    )

    theta = q0[1]

    assert_allclose(geod.trajectory[1][:, 1], theta, atol=1e-4, rtol=1e-4)


def test_schwarzschild_precession():
    q0 = [40., np.pi / 2, 0.]
    p0 = [0., 0., 3.83405]
    a = 0.
    end_lambda = 2000.
    step_size = 0.5

    geod = Timelike(
        position=q0,
        momentum=p0,
        a=a,
        end_lambda=end_lambda,
        step_size=step_size,
        return_cartesian=False,
        julia=True
    )

    assert_allclose(geod.trajectory[1][:, 1], np.pi/2, atol=1e-12, rtol=1e-12)
    assert_allclose(geod.trajectory[1][:, 4], 0, atol=1e-12, rtol=1e-12)


def test_kerr_frame_dragging():
    """
    Tests, if higher spin implies a "faster" capture (in terms of lambda),
    owing to frame dragging effects

    """
    q0 = [2.5, np.pi / 2, 0.]
    p0 = [0., 0., -8.5]
    end_lambda = 10.
    step_size = 0.005

    sch_geod = Timelike(
        position=q0,
        momentum=p0,
        a=0.,
        end_lambda=end_lambda,
        step_size=step_size,
        return_cartesian=False,
        julia=False
    )

    kerr_geod = Timelike(
        position=q0,
        momentum=p0,
        a=0.9,
        end_lambda=end_lambda,
        step_size=step_size,
        return_cartesian=False,
        julia=False
    )

    sch_traj = sch_geod.trajectory
    kerr_traj = kerr_geod.trajectory

    # Final lambda_sch > Final lambda_kerr
    assert sch_traj[0].shape[0] > kerr_traj[0].shape[0]
    assert sch_traj[0][-1] > kerr_traj[0][-1]
    # Final r_sch > Final r_kerr
    assert sch_traj[1][:, 0][-1] > kerr_traj[1][:, 0][-1]
