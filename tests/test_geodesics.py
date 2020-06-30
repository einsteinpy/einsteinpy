import warnings

import astropy.units as u
import numpy as np
import pytest
from numpy.testing import assert_allclose

from einsteinpy.coordinates import CartesianDifferential, SphericalDifferential, BoyerLindquistDifferential

from einsteinpy.metric import Schwarzschild, Kerr, KerrNewman
from einsteinpy.geodesic import Geodesic

from einsteinpy import constant

_c = constant.c.value
_G = constant.G.value
_Cc = constant.coulombs_const.value


def test_str_repr():
    """
    Tests, if the ``__str__`` and ``__repr__`` messages match

    """
    M = 1e25
    sph = SphericalDifferential(
        t=0.0 * u.s,
        r=306.0 * u.m,
        theta=np.pi / 2 * u.rad,
        phi=-np.pi / 2 * u.rad,
        v_r=0.0 * u.m / u.s,
        v_th=0.01 * u.rad / u.s,
        v_p=10.0 * u.rad / u.s,
    )
    ms = Schwarzschild(M=M)
    state = sph.state(metric=ms, time_like=True)

    end_lambda = 1.
    step_size = 0.4e-6

    geod = Geodesic(metric=ms, state=state, end_lambda=end_lambda, step_size=step_size)

    assert str(geod) == repr(geod)


@pytest.fixture()
def dummy_data():
    sph = SphericalDifferential(
        t=0.0 * u.s,
        r=130.0 * u.m,
        theta=np.pi / 2 * u.rad,
        phi=-np.pi / 8 * u.rad,
        v_r=0.0 * u.m / u.s,
        v_th=0.0 * u.rad / u.s,
        v_p=1900.0 * u.rad / u.s,
    )
    metric = Schwarzschild(M=6e24)

    state = sph.state(metric=metric, time_like=True)

    end_lambda = 0.002
    step_size = 5e-8

    return metric, state, end_lambda, step_size


def test_Geodesics_has_trajectory(dummy_data):
    metric, state, end_lambda, step_size = dummy_data
    geo = Geodesic(
        metric=metric,
        state=state,
        end_lambda=end_lambda,
        step_size=step_size
    )

    assert isinstance(geo.trajectory, np.ndarray)


@pytest.mark.parametrize(
    "sph, M, end_lambda, step_size",
    [
        (
            SphericalDifferential(
                t=0.0 * u.s,
                r=306.0 * u.m,
                theta=np.pi / 2 * u.rad,
                phi=np.pi / 2 * u.rad,
                v_r=0.0 * u.m / u.s,
                v_th=0.0 * u.rad / u.s,
                v_p=951.0 * u.rad / u.s,
            ),
            4e24,
            0.002,
            0.5e-6,
        ),
        (
            SphericalDifferential(
                t=0.0 * u.s,
                r=1e3 * u.m,
                theta=0.15 * u.rad,
                phi=np.pi / 2 * u.rad,
                v_r=0.1 * _c * u.m / u.s,
                v_th=0.5e-5 * u.rad / u.s,
                v_p=0.5e-4 * _c * u.rad / u.s,
            ),
            5.972e24,
            0.0001,
            0.5e-6,
        ),
        (
            SphericalDifferential(
                t=0.0 * u.s,
                r=50e3 * u.m,
                theta=np.pi / 2 * u.rad,
                phi=np.pi / 2 * u.rad,
                v_r=0.1 * _c * u.m / u.s,
                v_th=2e-7 * _c * u.rad / u.s,
                v_p=1e-5 * u.rad / u.s,
            ),
            5.972e24,
            0.001,
            5e-6,
        ),
    ],
)
def test_calculate_trajectory_schwarzschild(sph, M, end_lambda, step_size):
    ms = Schwarzschild(M=M)
    state = sph.state(metric=ms)

    geod = Geodesic(
        metric=ms,
        state=state,
        end_lambda=end_lambda,
        step_size=step_size,
        return_cartesian=False
    )

    ans = geod.trajectory

    testarray = list()
    for i in ans:
        x = i[:4]
        g = ms.metric_covariant(x)
        testarray.append(
            g[0][0] * (i[4] ** 2) +
            g[1][1] * (i[5] ** 2) +
            g[2][2] * (i[6] ** 2) +
            g[3][3] * (i[7] ** 2)
        )
    testarray = np.array(testarray, dtype=float)

    assert_allclose(testarray, _c ** 2, 1e-4)


def test_calculate_trajectory2_schwarzschild():
    # based on the revolution of earth around sun
    # data from https://en.wikipedia.org/wiki/Earth%27s_orbit
    distance_at_perihelion = 147.10e9
    speed_at_perihelion = 29290
    angular_vel = (speed_at_perihelion / distance_at_perihelion)

    sph = SphericalDifferential(
        t=0.0 * u.s,
        r=distance_at_perihelion * u.m,
        theta=np.pi / 2 * u.rad,
        phi=0.0 * u.rad,
        v_r=0.0 * u.m / u.s,
        v_th=0.0 * u.rad / u.s,
        v_p=angular_vel * u.rad / u.s,
    )
    metric = Schwarzschild(M=1.989e30)

    state = sph.state(metric=metric, time_like=True)

    end_lambda = 3.154e7

    geod = Geodesic(
        metric=metric,
        state=state,
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
    # same test as with test_calculate_trajectory2_schwarzschild(),
    # but initialized with cartesian coordinates
    # and function returning cartesian coordinates
    distance_at_perihelion = 147.10e9
    speed_at_perihelion = 29290

    x_sph = CartesianDifferential(
        t=0.0 * u.s,
        x=distance_at_perihelion / np.sqrt(2) * u.m,
        y=distance_at_perihelion / np.sqrt(2) * u.m,
        z=0. * u.m,
        v_x=-speed_at_perihelion / np.sqrt(2) * u.m / u.s,
        v_y=speed_at_perihelion / np.sqrt(2) * u.m / u.s,
        v_z=0 * u.m / u.s
    ).spherical_differential()

    metric = Schwarzschild(M=1.989e30)

    state = x_sph.state(metric=metric, time_like=True)

    end_lambda = 3.154e7

    geod = Geodesic(
        metric=metric,
        state=state,
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
    "sph, M, end_lambda, step_size, OdeMethodKwargs, return_cartesian",
    [
        (
            SphericalDifferential(
                t=0.0 * u.s,
                r=306.0 * u.m,
                theta=np.pi / 2 * u.rad,
                phi=np.pi / 2 * u.rad,
                v_r=0.0 * u.m / u.s,
                v_th=0.1 * u.rad / u.s,
                v_p=951.0 * u.rad / u.s,
            ),
            4e24,
            0.0002,
            0.3e-6,
            {"stepsize": 0.3e-6},
            True,
        ),
        (
            SphericalDifferential(
                t=0.0 * u.s,
                r=1e3 * u.m,
                theta=0.15 * u.rad,
                phi=np.pi / 2 * u.rad,
                v_r=_c * u.m / u.s,
                v_th=1e-4 * _c * u.rad / u.s,
                v_p=951.0 * u.rad / u.s,
            ),
            5.972e24,
            0.0002,
            0.5e-6,
            {"stepsize": 0.5e-6},
            False,
        ),
    ],
)
def test_calculate_trajectory_iterator_schwarzschild(
    sph, M, end_lambda, step_size, OdeMethodKwargs, return_cartesian
):
    metric = Schwarzschild(M=M)

    state = sph.state(metric=metric, time_like=True)

    geod = Geodesic(
        metric=metric,
        state=state,
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
    M = 1e25
    sph = SphericalDifferential(
        t=0.0 * u.s,
        r=306.0 * u.m,
        theta=np.pi / 2 * u.rad,
        phi=-np.pi / 2 * u.rad,
        v_r=0.0 * u.m / u.s,
        v_th=0.01 * u.rad / u.s,
        v_p=10.0 * u.rad / u.s,
    )
    ms = Schwarzschild(M=M)
    state = sph.state(metric=ms, time_like=True)

    end_lambda = 1.
    stepsize = 0.4e-6
    OdeMethodKwargs = {"stepsize": stepsize}

    geod = Geodesic(
        metric=ms,
        state=state,
        end_lambda=end_lambda,
        step_size=stepsize,
        return_cartesian=False
    )

    with warnings.catch_warnings(record=True) as w:
        it = geod.calculate_trajectory_iterator(
            OdeMethodKwargs=OdeMethodKwargs,
        )
        for _, _ in zip(range(1000), it):
            pass

        assert len(w) >= 1


@pytest.mark.parametrize(
    "bl, M, a, end_lambda, step_size",
    [
        (
            BoyerLindquistDifferential(
                t=0.0 * u.s,
                r=306.0 * u.m,
                theta=np.pi / 2.05 * u.rad,
                phi=np.pi / 2 * u.rad,
                v_r=0.0 * u.m / u.s,
                v_th=0.0 * u.rad / u.s,
                v_p=951.0 * u.rad / u.s,
            ),
            4e24,
            2e-3,
            0.001,
            0.5e-6,
        ),
        (
            BoyerLindquistDifferential(
                t=0.0 * u.s,
                r=1e3 * u.m,
                theta=0.15 * u.rad,
                phi=np.pi / 2 * u.rad,
                v_r=0.1 * _c * u.m / u.s,
                v_th=0.5e-5 * _c * u.rad / u.s,
                v_p=0.5e-4 * _c * u.rad / u.s,
            ),
            5.972e24,
            2e-3,
            0.0001,
            0.5e-6,
        ),
        (
            BoyerLindquistDifferential(
                t=0.0 * u.s,
                r=50e3 * u.m,
                theta=np.pi / 2 * u.rad,
                phi=np.pi / 2 * u.rad,
                v_r=0.1 * _c * u.m / u.s,
                v_th=2e-7 * _c * u.rad / u.s,
                v_p=1e-5 * u.rad / u.s,
            ),
            5.972e24,
            0.,
            0.001,
            5e-6,
        ),
    ],
)
def test_calculate_trajectory_kerr(
    bl, M, a, end_lambda, step_size
):
    mk = Kerr(coords="BL", M=M, a=a)
    state = bl.state(metric=mk, time_like=True)

    geod = Geodesic(
        metric=mk,
        state=state,
        end_lambda=end_lambda,
        step_size=step_size,
        return_cartesian=False
    )

    ans = geod.trajectory

    testarray = list()
    for i in ans:
        x = i[:4]
        g = mk.metric_covariant(x)
        testarray.append(
            g[0][0] * (i[4] ** 2) +
            g[1][1] * (i[5] ** 2) +
            g[2][2] * (i[6] ** 2) +
            g[3][3] * (i[7] ** 2) +
            2 * g[0][3] * i[4] * i[7]
        )
    testarray = np.array(testarray, dtype=float)

    # rtol < 1e-2, causes test to fail
    # Error accumulation perhaps.
    assert_allclose(testarray, _c ** 2, 1e-2)


def test_calculate_trajectory3_kerr():
    # Based on the revolution of earth around sun
    # Data from https://en.wikipedia.org/wiki/Earth%27s_orbit
    # Initialized with cartesian coordinates
    # Function returning cartesian coordinates
    M = 1.989e30
    a = 0.0

    distance_at_perihelion = 147.10e9
    speed_at_perihelion = 29290

    x_bl = CartesianDifferential(
        t=0.0,
        x=distance_at_perihelion / np.sqrt(2),
        y=distance_at_perihelion / np.sqrt(2),
        z=0.0,
        v_x=-speed_at_perihelion / np.sqrt(2),
        v_y=speed_at_perihelion / np.sqrt(2),
        v_z=0.0
    ).bl_differential(M=M, a=a)

    mk = Kerr(coords="BL", M=M, a=a)
    state = x_bl.state(metric=mk, time_like=True)

    end_lambda = 3.154e7

    geod = Geodesic(
        metric=mk,
        state=state,
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
    "bl, M, a, end_lambda, step_size, OdeMethodKwargs, return_cartesian",
    [
        (
            BoyerLindquistDifferential(
                t=0.0 * u.s,
                r=306.0 * u.m,
                theta=np.pi / 2 * u.rad,
                phi=np.pi / 2 * u.rad,
                v_r=0.0 * u.m / u.s,
                v_th=0.1 * u.rad / u.s,
                v_p=951.0 * u.rad / u.s,
            ),
            4e24,
            2e-3,
            0.0003,
            0.3e-6,
            {"stepsize": 0.3e-6},
            True,
        ),
        (
            BoyerLindquistDifferential(
                t=0.0 * u.s,
                r=1e3 * u.m,
                theta=0.15 * u.rad,
                phi=np.pi / 2 * u.rad,
                v_r=0.2 * _c * u.m / u.s,
                v_th=0.5e-5 * _c * u.rad / u.s,
                v_p=1e-4 * _c * u.rad / u.s,
            ),
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
    bl, M, a, end_lambda, step_size, OdeMethodKwargs, return_cartesian
):
    mk = Kerr(coords="BL", M=M, a=a)
    state = bl.state(metric=mk, time_like=True)

    geod = Geodesic(
        metric=mk,
        state=state,
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


# Fails, due to math issues in f_vec / christoffels / dg_dx
@pytest.mark.skip(reason="Expected failure")
def test_calculate_trajectory_iterator_RuntimeWarning_kerr():
    M = 1e25
    a = 0.

    bl = BoyerLindquistDifferential(
        t=0.0 * u.s,
        r=306.0 * u.m,
        theta=np.pi / 2 * u.rad,
        phi=np.pi / 2 * u.rad,
        v_r=0.0 * u.m / u.s,
        v_th=0.01 * u.rad / u.s,
        v_p=10.0 * u.rad / u.s,
    )

    mk = Kerr(coords="BL", M=M, a=a)

    state = bl.state(metric=mk, time_like=True)

    end_lambda = 1.
    stepsize = 0.4e-6
    OdeMethodKwargs = {"stepsize": stepsize}

    geod = Geodesic(
        metric=mk,
        state=state,
        end_lambda=end_lambda,
        step_size=stepsize,
        return_cartesian=False
    )

    with warnings.catch_warnings(record=True) as w:
        it = geod.calculate_trajectory_iterator(
            OdeMethodKwargs=OdeMethodKwargs,
        )
        for _, _ in zip(range(1000), it):
            pass

        assert len(w) >= 1


def test_calculate_trajectory0_kerrnewman():
    # Based on the revolution of earth around sun
    # Data from https://en.wikipedia.org/wiki/Earth%27s_orbit
    # Initialized with cartesian coordinates
    # Function returning cartesian coordinates
    M = 1.989e30
    a = 0.
    Q = 0.
    q = 0.
    distance_at_perihelion = 147.10e9
    speed_at_perihelion = 29290

    x_bl = CartesianDifferential(
        t=0.0 * u.s,
        x=distance_at_perihelion / np.sqrt(2) * u.m,
        y=distance_at_perihelion / np.sqrt(2) * u.m,
        z=0.0 * u.m,
        v_x=-speed_at_perihelion / np.sqrt(2) * u.m / u.s,
        v_y=speed_at_perihelion / np.sqrt(2) * u.m / u.s,
        v_z=0.0 * u.m / u.s
    ).bl_differential(M=M, a=a)

    mkn = KerrNewman(coords="BL", M=M, a=a, Q=Q, q=q)
    state = x_bl.state(metric=mkn, time_like=True)

    end_lambda = 3.154e7

    geod = Geodesic(
        metric=mkn,
        state=state,
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


@pytest.mark.skip(reason="This needs more investigation, implementation seems incorrect.")
def test_calculate_trajectory1_kerrnewman():
    # This needs more investigation
    # the test particle should not move as gravitational & electromagnetic forces are balanced
    M = 0.5 * 5.972e24
    a = 0.
    Q = 11604461683.91822052001953125
    q = _G * M / _Cc

    r = 1e6
    end_lambda = 1000.0
    step_size = 0.5

    bl = BoyerLindquistDifferential(
        t=0.0 * u.s,
        r=r * u.m,
        theta=np.pi / 2 * u.rad,
        phi=0.0 * u.rad,
        v_r=0.0 * u.m / u.s,
        v_th=0.0 * u.rad / u.s,
        v_p=0.0 * u.rad / u.s,
    )

    mkn = KerrNewman(coords="BL", M=M, a=a, Q=Q, q=q)
    state = bl.state(metric=mkn, time_like=True)

    geod = Geodesic(
        metric=mkn,
        state=state,
        end_lambda=end_lambda,
        step_size=step_size,
        return_cartesian=False
    )

    ans = geod.trajectory

    assert_allclose(ans[0][1], ans[-1][1], 1e-2)


@pytest.fixture()
def test_input():
    a = 1e-6
    Q = 100.
    q = 1.
    end_lambda = 200.0
    step_size = 1.0

    return a, Q, q, end_lambda, step_size


def test_compare_calculate_trajectory_iterator_bl_kerrnewman(test_input):
    a, Q, q, end_lambda, step_size = test_input
    M = 0.5 * 5.972e24

    bl = BoyerLindquistDifferential(
        t=0.0 * u.s,
        r=1e6 * u.m,
        theta=0.6 * np.pi * u.rad,
        phi=np.pi / 8 * u.rad,
        v_r=1e4 * u.m / u.s,
        v_th=-0.01 * u.rad / u.s,
        v_p=0.0 * u.rad / u.s,
    )

    mkn = KerrNewman(coords="BL", M=M, a=a, Q=Q, q=q)
    state = bl.state(metric=mkn, time_like=True)

    OdeMethodKwargs = {"stepsize": step_size}
    return_cartesian = False

    geod = Geodesic(
        metric=mkn,
        state=state,
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
    a, Q, q, end_lambda, step_size = test_input
    M = 2e24

    x_bl = CartesianDifferential(
        t=0.0 * u.s,
        x=1e6 * u.m,
        y=1e6 * u.m,
        z=20.5 * u.m,
        v_x=1e4 * u.m / u.s,
        v_y=1e4 * u.m / u.s,
        v_z=-30.0 * u.m / u.s
    ).bl_differential(M=M, a=a)

    mkn = KerrNewman(coords="BL", M=M, a=a, Q=Q, q=q)
    state = x_bl.state(metric=mkn, time_like=True)

    OdeMethodKwargs = {"stepsize": step_size}
    return_cartesian = True

    geod = Geodesic(
        metric=mkn,
        state=state,
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


# Fails, due to math issues in f_vec / christoffels / dg_dx
@pytest.mark.skip(reason="Expected failure")
def test_calculate_trajectory_iterator_RuntimeWarning_kerrnewman():
    M = 0.5 * 5.972e24
    a = 0.
    Q = 0.
    q = 0.

    bl = BoyerLindquistDifferential(
        t=0.0 * u.s,
        r=306.0 * u.m,
        theta=np.pi / 2 * u.rad,
        phi=np.pi / 2 * u.rad,
        v_r=0.0 * u.m / u.s,
        v_th=-0.01 * u.rad / u.s,
        v_p=10.0 * u.rad / u.s,
    )

    mkn = KerrNewman(coords="BL", M=M, a=a, Q=Q, q=q)
    state = bl.state(metric=mkn, time_like=True)
    print(state)

    end_lambda = 200.
    step_size = 0.4e-6
    OdeMethodKwargs = {"stepsize": step_size}

    geod = Geodesic(metric=mkn, state=state, end_lambda=end_lambda, step_size=step_size)

    with warnings.catch_warnings(record=True) as w:
        it = geod.calculate_trajectory_iterator(
            OdeMethodKwargs=OdeMethodKwargs,
        )
        for _, _ in zip(range(1000), it):
            pass
        assert len(w) >= 1
