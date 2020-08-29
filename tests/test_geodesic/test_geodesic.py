import warnings

import numpy as np
import pytest
from numpy.testing import assert_allclose
from unittest.mock import Mock

from einsteinpy.examples import precession
from einsteinpy.geodesic import Geodesic, Nulllike, Timelike


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


def test_str_repr(dummy_time_python):
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


def test_geodesic_attribute1(dummy_time_python):
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
    traj = geod.trajectory

    assert isinstance(traj, tuple)
    assert isinstance(traj[0], np.ndarray)
    assert isinstance(traj[1], np.ndarray)


def test_geodesic_attribute2(dummy_time_python):
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
    traj = geod.trajectory

    assert traj
    assert traj[0].shape[0] == traj[1].shape[0]
    assert traj[1].shape[1] == 6


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


partial_lambdas = np.array([0.0000e+00, 5.0000e-01, 1.0000e+00, 1.9990e+03, 1.9995e+03, 2.0000e+03])
partial_vecs = np.array(
    [[ 4.00000000e+01,  0.00000000e+00,  2.44929360e-15, 0.00000000e+00,  0.00000000e+00,  3.83405000e+00],
    [ 3.99999197e+01, -4.79255517e-02,  2.44929044e-15, 2.17126283e-04, -2.81284902e-19,  3.83405000e+00],
    [ 3.99996789e+01, -9.58508493e-02,  2.44928097e-15, 4.34253495e-04, -5.62571255e-19,  3.83405000e+00],
    [-6.00084085e+00, -2.52448871e+00,  3.98636802e-16, -2.33673001e-01, -9.49738949e-15,  3.83405000e+00],
    [-6.18465627e+00, -2.28046143e+00,  4.03625070e-16, -2.36036036e-01, -9.50787744e-15,  3.83405000e+00],
    [-6.35745272e+00, -2.03236148e+00,  4.08689510e-16, -2.38256302e-01, -9.51810736e-15,  3.83405000e+00]]
)

precession = Mock()
precession.return_value = (partial_lambdas, partial_vecs)

def test_precession_attr():
    """
    Checks various attributes of precession()

    """
    p = precession()
    L = p[1][0, 5]

    assert p[0].shape[0] == p[1].shape[1]
    assert p[1].shape[1] == 6
    assert_allclose(p[1][:, 4], 0, atol=1e-12, rtol=1e-12)
    assert_allclose(p[1][:, 5], L, atol=1e-12, rtol=1e-12)
