import warnings

import astropy.units as u
import numpy as np
import pytest
from numpy.testing import assert_allclose

from einsteinpy.metric import Schwarzschild, Kerr, KerrNewman
from einsteinpy.coordinates import CartesianConversion
from einsteinpy.coordinates.utils import four_position, stacked_vec
from einsteinpy.geodesic import Geodesic

from einsteinpy import constant

_c = constant.c.value


@pytest.fixture()
def dummy_data():
    M = 6e24
    t = 0.
    x_vec = np.array([130.0, np.pi / 2, -np.pi / 8])
    v_vec = np.array([0.0, 0.0, 1900.0])

    metric = Schwarzschild(M=M)
    x_4vec = four_position(t, x_vec)
    metric_mat = metric.metric_covariant(x_4vec)
    init_vec = stacked_vec(metric_mat, t, x_vec, v_vec, time_like=True)

    end_lambda = 0.002
    step_size = 5e-8

    return metric, init_vec, end_lambda, step_size


def test_Geodesics_has_trajectory(dummy_data):
    metric, init_vec, end_lambda, step_size = dummy_data
    geo = Geodesic(
        metric=metric,
        init_vec=init_vec,
        end_lambda=end_lambda,
        step_size=step_size
    )

    assert isinstance(geo.trajectory, np.ndarray)


@pytest.mark.parametrize(
    "x_vec, v_vec, t, M, end_lambda, step_size",
    [
        (
            np.array([306., np.pi / 2, np.pi / 2]),
            np.array([0., 0., 951.]),
            0.,
            4e24,
            0.002,
            0.5e-6,
        ),
        (
            np.array([1e3, 0.15, np.pi / 2]),
            np.array([0.1 * _c, 0.5e-5 * _c, 0.5e-4 * _c]),
            0.,
            5.972e24,
            0.0001,
            0.5e-6,
        ),
        (
            np.array([50e3, np.pi / 2, np.pi / 2]),
            np.array([0.1 * _c, 2e-7 * _c, 1e-5]),
            0.,
            5.972e24,
            0.001,
            5e-6,
        ),
    ],
)
def test_calculate_trajectory_schwarzschild(
    x_vec, v_vec, t, M, end_lambda, step_size
):
    ms_cov = Schwarzschild(M=M)
    x_4vec = four_position(t, x_vec)
    ms_cov_mat = ms_cov.metric_covariant(x_4vec)
    init_vec = stacked_vec(ms_cov_mat, t, x_vec, v_vec, time_like=True)

    geod = Geodesic(
        metric=ms_cov,
        init_vec=init_vec,
        end_lambda=end_lambda,
        step_size=step_size,
        return_cartesian=False
    )

    ans = geod.trajectory

    testarray = list()
    for i in ans:
        x = i[:4]
        g = ms_cov.metric_covariant(x)
        testarray.append(
            g[0][0] * (i[4] ** 2) +
            g[1][1] * (i[5] ** 2) +
            g[2][2] * (i[6] ** 2) +
            g[3][3] * (i[7] ** 2)
        )
    testarray = np.array(testarray, dtype=float)

    assert_allclose(testarray, 1., 1e-4)


def test_calculate_trajectory2_schwarzschild():
    # based on the revolution of earth around sun
    # data from https://en.wikipedia.org/wiki/Earth%27s_orbit
    t = 0.
    M = 1.989e30
    distance_at_perihelion = 147.10e9
    speed_at_perihelion = 30290
    angular_vel = (speed_at_perihelion / distance_at_perihelion)

    x_vec = np.array([distance_at_perihelion, np.pi / 2, 0])
    v_vec = np.array([0.0, 0.0, angular_vel])

    ms_cov = Schwarzschild(M=M)
    x_4vec = four_position(t, x_vec)
    ms_cov_mat = ms_cov.metric_covariant(x_4vec)
    init_vec = stacked_vec(ms_cov_mat, t, x_vec, v_vec, time_like=True)

    end_lambda = 3.154e7

    geod = Geodesic(
        metric=ms_cov,
        init_vec=init_vec,
        end_lambda=end_lambda,
        step_size=end_lambda / 2e3,
        return_cartesian=False
    )

    ans = geod.trajectory

    # velocity should be 29.29 km/s at aphelion(where r is max)
    i = np.argmax(ans[:, 1])  # index where radial distance is max
    v_aphelion = (((ans[i][1] * ans[i][7]) * (u.m / u.s)).to(u.km / u.s)).value

    assert_allclose(v_aphelion, 29.29, rtol=0.01)


def test_calculate_trajectory3_schwarzschild():
    # same test as with test_calculate_trajectory2(),
    # but initialized with cartesian coordinates
    # and function returning cartesian coordinates
    t = 0.
    M = 1.989e30
    distance_at_perihelion = 147.10e9
    speed_at_perihelion = 30290

    x_sph = CartesianConversion(distance_at_perihelion / np.sqrt(2), distance_at_perihelion / np.sqrt(2), 0., -speed_at_perihelion / np.sqrt(2), speed_at_perihelion / np.sqrt(2), 0.).convert_spherical()

    x_vec = x_sph[:3]
    v_vec = x_sph[3:]

    ms_cov = Schwarzschild(M=M)
    x_4vec = four_position(t, x_vec)
    ms_cov_mat = ms_cov.metric_covariant(x_4vec)
    init_vec = stacked_vec(ms_cov_mat, t, x_vec, v_vec, time_like=True)

    end_lambda = 3.154e7

    geod = Geodesic(
        metric=ms_cov,
        init_vec=init_vec,
        end_lambda=end_lambda,
        step_size=end_lambda / 2e3,
    )

    ans = geod.trajectory

    # velocity should be 29.29 km/s at aphelion(where r is max)
    R = np.sqrt(ans[:, 1] ** 2 + ans[:, 2] ** 2 + ans[:, 3] ** 2)
    i = np.argmax(R)  # index where radial distance is max
    v_aphelion = (
        (np.sqrt(ans[i, 5] ** 2 + ans[i, 6] ** 2 + ans[i, 7] ** 2) * (u.m / u.s)).to(
            u.km / u.s
        )
    ).value

    assert_allclose(v_aphelion, 29.29, rtol=0.01)


@pytest.mark.parametrize(
    "x_vec, v_vec, t, M, end_lambda, step_size, OdeMethodKwargs, return_cartesian",
    [
        (
            np.array([306., np.pi / 2, np.pi / 2]),
            np.array([0., 0.1, 951.]),
            0.,
            4e24,
            0.0002,
            0.3e-6,
            {"stepsize": 0.3e-6},
            True,
        ),
        (
            np.array([1e3, 0.15, np.pi / 2]),
            np.array([_c, 0.5e-5 * _c, 1e-4 * _c]),
            0.,
            5.972e24,
            0.0002,
            0.5e-6,
            {"stepsize": 0.5e-6},
            False,
        ),
    ],
)
def test_calculate_trajectory_iterator_schwarzschild(
    x_vec, v_vec, t, M, end_lambda, step_size, OdeMethodKwargs, return_cartesian
):
    ms_cov = Schwarzschild(M=M)
    x_4vec = four_position(t, x_vec)
    ms_cov_mat = ms_cov.metric_covariant(x_4vec)
    init_vec = stacked_vec(ms_cov_mat, t, x_vec, v_vec, time_like=True)

    geod = Geodesic(
        metric=ms_cov,
        init_vec=init_vec,
        end_lambda=end_lambda,
        step_size=step_size,
        return_cartesian=return_cartesian
    )

    traj = geod.trajectory

    traj_iter = geod.calculate_trajectory_iterator(OdeMethodKwargs=OdeMethodKwargs, return_cartesian=return_cartesian)

    traj_iter_list = list()
    for _, val in zip(range(50), traj_iter):
        traj_iter_list.append(val[1])
    traj_iter_arr = np.array(traj_iter_list)

    assert_allclose(traj[:50, :], traj_iter_arr, rtol=1e-10)


def test_calculate_trajectory_iterator_RuntimeWarning_schwarzschild():
    t = 0.
    M = 1e25

    x_vec = np.array([306., np.pi / 2, np.pi / 2])
    v_vec = np.array([0., 0.01, 10.])

    ms_cov = Schwarzschild(M=M)
    x_4vec = four_position(t, x_vec)
    ms_cov_mat = ms_cov.metric_covariant(x_4vec)
    init_vec = stacked_vec(ms_cov_mat, t, x_vec, v_vec, time_like=True)

    end_lambda = 1.
    stepsize = 0.4e-6
    OdeMethodKwargs = {"stepsize": stepsize}

    geod = Geodesic(
        metric=ms_cov,
        init_vec=init_vec,
        end_lambda=end_lambda,
        step_size=stepsize,
        return_cartesian=False
    )

    with warnings.catch_warnings(record=True) as w:
        it = geod.calculate_trajectory_iterator(
            OdeMethodKwargs=OdeMethodKwargs,
            stop_on_singularity=True,
        )
        for _, _ in zip(range(1000), it):
            pass

        assert len(w) >= 1


def test_calculate_trajectory_iterator_RuntimeWarning2_schwarzschild():
    t = 0.
    M = 1e25

    x_vec = np.array([306., np.pi / 2, np.pi / 3])
    v_vec = np.array([0., 0.01, 10.])

    ms_cov = Schwarzschild(M=M)
    x_4vec = four_position(t, x_vec)
    ms_cov_mat = ms_cov.metric_covariant(x_4vec)
    init_vec = stacked_vec(ms_cov_mat, t, x_vec, v_vec, time_like=True)

    end_lambda = 1.
    stepsize = 0.4e-6
    OdeMethodKwargs = {"stepsize": stepsize}

    geod = Geodesic(
        metric=ms_cov,
        init_vec=init_vec,
        end_lambda=end_lambda,
        step_size=stepsize,
        return_cartesian=False
    )

    with warnings.catch_warnings(record=True) as w:
        it = geod.calculate_trajectory_iterator(
            OdeMethodKwargs=OdeMethodKwargs,
            stop_on_singularity=False,
        )
        for _, _ in zip(range(1000), it):
            pass

        assert len(w) >= 1


@pytest.mark.parametrize(
    "x_vec, v_vec, t, M, a, end_lambda, step_size",
    [
        (
            np.array([306., np.pi / 2.05, np.pi / 2]),
            np.array([0., 0., 951.]),
            0.,
            4e24,
            2e-3,
            0.001,
            0.5e-6,
        ),
        (
            np.array([1e3, 0.15, np.pi / 2]),
            np.array([0.1 * _c, 0.5e-5 * _c, 0.5e-4 * _c]),
            0.,
            5.972e24,
            2e-3,
            0.0001,
            0.5e-6,
        ),
        (
            np.array([50e3, np.pi / 2, np.pi / 2]),
            np.array([0.1 * _c, 2e-7 * _c, 1e-5]),
            0.,
            5.972e24,
            0.,
            0.001,
            5e-6,
        ),
    ],
)
def test_calculate_trajectory_kerr(
    x_vec, v_vec, t, M, a, end_lambda, step_size
):
    mk_cov = Kerr(coords="BL", M=M, a=a)
    x_4vec = four_position(t, x_vec)
    mk_cov_mat = mk_cov.metric_covariant(x_4vec)
    init_vec = stacked_vec(mk_cov_mat, t, x_vec, v_vec, time_like=True)

    geod = Geodesic(
        metric=mk_cov,
        init_vec=init_vec,
        end_lambda=end_lambda,
        step_size=step_size,
        return_cartesian=False
    )

    ans = geod.trajectory

    testarray = list()
    for i in ans:
        x = i[:4]
        g = mk_cov.metric_covariant(x)
        testarray.append(
            g[0][0] * (i[4] ** 2) +
            g[1][1] * (i[5] ** 2) +
            g[2][2] * (i[6] ** 2) +
            g[3][3] * (i[7] ** 2) +
            2 * g[0][3] * i[4] * i[7]
        )
    testarray = np.array(testarray, dtype=float)

    assert_allclose(testarray, 1., 1e-4)


def test_calculate_trajectory3_kerr():
    # Based on the revolution of earth around sun
    # Data from https://en.wikipedia.org/wiki/Earth%27s_orbit
    # Initialized with cartesian coordinates
    # Function returning cartesian coordinates
    t = 0.
    M = 1.989e30
    a = 0.
    distance_at_perihelion = 147.10e9
    speed_at_perihelion = 30290

    x_sph = CartesianConversion(distance_at_perihelion / np.sqrt(2), distance_at_perihelion / np.sqrt(2), 0., -speed_at_perihelion / np.sqrt(2), speed_at_perihelion / np.sqrt(2), 0.).convert_spherical()

    x_vec = x_sph[:3]
    v_vec = x_sph[3:]

    mk_cov = Kerr(coords="BL", M=M, a=a)
    x_4vec = four_position(t, x_vec)
    mk_cov_mat = mk_cov.metric_covariant(x_4vec)
    init_vec = stacked_vec(mk_cov_mat, t, x_vec, v_vec, time_like=True)

    end_lambda = 3.154e7

    geod = Geodesic(
        metric=mk_cov,
        init_vec=init_vec,
        end_lambda=end_lambda,
        step_size=end_lambda / 2e3,
    )

    ans = geod.trajectory

    # velocity should be 29.29 km/s at aphelion(where r is max)
    R = np.sqrt(ans[:, 1] ** 2 + ans[:, 2] ** 2 + ans[:, 3] ** 2)
    i = np.argmax(R)  # index where radial distance is max
    v_aphelion = (
        (np.sqrt(ans[i, 5] ** 2 + ans[i, 6] ** 2 + ans[i, 7] ** 2) * (u.m / u.s)).to(
            u.km / u.s
        )
    ).value

    assert_allclose(v_aphelion, 29.29, rtol=0.01)


@pytest.mark.parametrize(
    "x_vec, v_vec, t, M, a, end_lambda, step_size, OdeMethodKwargs, return_cartesian",
    [
        (
            np.array([306., np.pi / 2, np.pi / 2]),
            np.array([0., 0.1, 951.]),
            0.,
            4e24,
            2e-3,
            0.0003,
            0.3e-6,
            {"stepsize": 0.3e-6},
            True,
        ),
        (
            np.array([1e3, 0.15, np.pi / 2]),
            np.array([0.2 * _c, 0.5e-5 * _c, 1e-4 * _c]),
            0.,
            5.972e24,
            0.,
            0.0004,
            0.5e-6,
            {"stepsize": 0.5e-6},
            False,
        ),
    ],
)
def test_calculate_trajectory_iterator_kerr(
    x_vec, v_vec, t, M, a, end_lambda, step_size, OdeMethodKwargs, return_cartesian
):

    mk_cov = Kerr(coords="BL", M=M, a=a)
    x_4vec = four_position(t, x_vec)
    mk_cov_mat = mk_cov.metric_covariant(x_4vec)
    init_vec = stacked_vec(mk_cov_mat, t, x_vec, v_vec, time_like=True)

    geod = Geodesic(
        metric=mk_cov,
        init_vec=init_vec,
        end_lambda=end_lambda,
        step_size=step_size,
        return_cartesian=return_cartesian
    )

    traj = geod.trajectory

    traj_iter = geod.calculate_trajectory_iterator(OdeMethodKwargs=OdeMethodKwargs, return_cartesian=return_cartesian)

    traj_iter_list = list()
    for _, val in zip(range(50), traj_iter):
        traj_iter_list.append(val[1])
    traj_iter_arr = np.array(traj_iter_list)

    assert_allclose(traj[:50, :], traj_iter_arr, rtol=1e-10)


def test_calculate_trajectory_iterator_RuntimeWarning_kerr():
    t = 0.
    M = 1e25
    a = 0.

    x_vec = np.array([306., np.pi / 2, np.pi / 2])
    v_vec = np.array([0., 0.01, 10.])

    mk_cov = Kerr(coords="BL", M=M, a=a)
    x_4vec = four_position(t, x_vec)
    mk_cov_mat = mk_cov.metric_covariant(x_4vec)
    init_vec = stacked_vec(mk_cov_mat, t, x_vec, v_vec, time_like=True)

    end_lambda = 1.
    stepsize = 0.4e-6
    OdeMethodKwargs = {"stepsize": stepsize}

    geod = Geodesic(
        metric=mk_cov,
        init_vec=init_vec,
        end_lambda=end_lambda,
        step_size=stepsize,
        return_cartesian=False
    )

    with warnings.catch_warnings(record=True) as w:
        it = geod.calculate_trajectory_iterator(
            OdeMethodKwargs=OdeMethodKwargs,
            stop_on_singularity=True,
        )
        for _, _ in zip(range(1000), it):
            pass

        assert len(w) >= 1


def test_calculate_trajectory0_kerrnewman():
    # Based on the revolution of earth around sun
    # Data from https://en.wikipedia.org/wiki/Earth%27s_orbit
    # Initialized with cartesian coordinates
    # Function returning cartesian coordinates
    t = 0.
    M = 1.989e30
    a = 0.
    Q = 0.
    q = 0.
    distance_at_perihelion = 147.10e9
    speed_at_perihelion = 30290

    x_sph = CartesianConversion(distance_at_perihelion / np.sqrt(2), distance_at_perihelion / np.sqrt(2), 0., -speed_at_perihelion / np.sqrt(2), speed_at_perihelion / np.sqrt(2), 0.).convert_spherical()

    x_vec = x_sph[:3]
    v_vec = x_sph[3:]

    mkn_cov = KerrNewman(coords="BL", M=M, a=a, Q=Q, q=q)
    x_4vec = four_position(t, x_vec)
    mkn_cov_mat = mkn_cov.metric_covariant(x_4vec)
    init_vec = stacked_vec(mkn_cov_mat, t, x_vec, v_vec, time_like=True)

    end_lambda = 3.154e7

    geod = Geodesic(
        metric=mkn_cov,
        init_vec=init_vec,
        end_lambda=end_lambda,
        step_size=end_lambda / 1.5e3,
    )

    ans = geod.trajectory

    # velocity should be 29.29 km/s at aphelion(where r is max)
    R = np.sqrt(ans[:, 1] ** 2 + ans[:, 2] ** 2 + ans[:, 3] ** 2)
    i = np.argmax(R)  # index where radial distance is max
    v_aphelion = (
        (np.sqrt(ans[i, 5] ** 2 + ans[i, 6] ** 2 + ans[i, 7] ** 2) * (u.m / u.s)).to(
            u.km / u.s
        )
    ).value

    assert_allclose(v_aphelion, 29.29, rtol=0.01)

def test_calculate_trajectory1_kerrnewman():
    # This needs more investigation
    # the test particle should not move as gravitational & electromagnetic forces are balanced
    t = 0.
    M = 0.5 * 5.972e24
    a = 0.
    Q = 11604461683.91822052001953125
    q = _G * M / _Cc

    r = 1e6
    end_lambda = 1000.0
    step_size = 0.5

    x_vec = np.array([r, 0.5 * np.pi, 0.])
    v_vec = np.array([0., 0., 0.])

    mkn_cov = KerrNewman(coords="BL", M=M, a=a, Q=Q, q=q)
    x_4vec = four_position(t, x_vec)
    mkn_cov_mat = mkn_cov.metric_covariant(x_4vec)
    init_vec = stacked_vec(mkn_cov_mat, t, x_vec, v_vec, time_like=True)

    geod = Geodesic(
        metric=mkn_cov,
        init_vec=init_vec,
        end_lambda=end_lambda,
        step_size=step_size,
        return_cartesian=False
    )

    ans = geod.trajectory

    assert_allclose(ans[0][1], ans[-1][1], 1e-2)


@pytest.fixture()
def test_input():
    t = 0.
    a = 1e-6
    Q = 100.
    q = 1.
    end_lambda = 200.0
    step_size = 1.0

    return t, a, Q, q, end_lambda, step_size


def test_compare_calculate_trajectory_iterator_bl_kerrnewman(test_input):
    t, a, Q, q, end_lambda, step_size = test_input
    M = 0.5 * 5.972e24

    x_vec = np.array([1e6, 0.6 * np.pi, np.pi / 8])
    v_vec = np.array([1e4, -0.01, 0.])

    mkn_cov = KerrNewman(coords="BL", M=M, a=a, Q=Q, q=q)
    x_4vec = four_position(t, x_vec)
    mkn_cov_mat = mkn_cov.metric_covariant(x_4vec)
    init_vec = stacked_vec(mkn_cov_mat, t, x_vec, v_vec, time_like=True)

    OdeMethodKwargs = {"stepsize": step_size}
    return_cartesian = False

    geod = Geodesic(
        metric=mkn_cov,
        init_vec=init_vec,
        end_lambda=end_lambda,
        step_size=step_size,
        return_cartesian=return_cartesian
    )

    traj = geod.trajectory

    traj_iter = geod.calculate_trajectory_iterator(OdeMethodKwargs=OdeMethodKwargs, return_cartesian=return_cartesian)

    traj_iter_list = list()
    for _, val in zip(range(20), traj_iter):
        traj_iter_list.append(val[1])
    traj_iter_arr = np.array(traj_iter_list)

    assert_allclose(traj[:20], traj_iter_arr)


def test_compare_calculate_trajectory_iterator_cartesian_kerrnewman(test_input):
    t, a, Q, q, end_lambda, step_size = test_input
    M = 2e24

    x_vec = np.array([1e6, 1e6, 20.5])
    v_vec = np.array([1e4, 1e4, -30.])

    mkn_cov = KerrNewman(coords="BL", M=M, a=a, Q=Q, q=q)
    x_4vec = four_position(t, x_vec)
    mkn_cov_mat = mkn_cov.metric_covariant(x_4vec)
    init_vec = stacked_vec(mkn_cov_mat, t, x_vec, v_vec, time_like=True)

    OdeMethodKwargs = {"stepsize": step_size}
    return_cartesian = True

    geod = Geodesic(
        metric=mkn_cov,
        init_vec=init_vec,
        end_lambda=end_lambda,
        step_size=step_size,
        return_cartesian=return_cartesian
    )

    traj = geod.trajectory

    traj_iter = geod.calculate_trajectory_iterator(OdeMethodKwargs=OdeMethodKwargs, return_cartesian=return_cartesian)

    traj_iter_list = list()
    for _, val in zip(range(20), traj_iter):
        traj_iter_list.append(val[1])
    traj_iter_arr = np.array(traj_iter_list)

    assert_allclose(traj[:20], traj_iter_arr)


def test_calculate_trajectory_iterator_RuntimeWarning_kerrnewman():
    M = 0.5 * 5.972e24
    a = 0.
    Q = 0.
    q = 0.
    t = 0.

    x_vec = np.array([306, np.pi / 2, np.pi / 2])
    v_vec = np.array([0., 0.01, 10.])

    mkn_cov = KerrNewman(coords="BL", M=M, a=a, Q=Q, q=q)
    x_4vec = four_position(t, x_vec)
    mkn_cov_mat = mkn_cov.metric_covariant(x_4vec)
    init_vec = stacked_vec(mkn_cov_mat, t, x_vec, v_vec, time_like=True)

    end_lambda = 200.
    step_size = 0.4e-6
    OdeMethodKwargs = {"stepsize": step_size}

    geod = Geodesic(metric=mkn_cov, init_vec=init_vec, end_lambda=end_lambda, step_size=step_size)

    with warnings.catch_warnings(record=True) as w:
        it = geod.calculate_trajectory_iterator(
            OdeMethodKwargs=OdeMethodKwargs,
            stop_on_singularity=True,
        )
        for _, _ in zip(range(1000), it):
            pass
        assert len(w) >= 1
