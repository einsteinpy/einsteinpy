import warnings

import astropy.units as u
import numpy as np
import pytest
from numpy.testing import assert_allclose

from einsteinpy.coordinates import CartesianDifferential, SphericalDifferential, BoyerLindquistDifferential

from einsteinpy.metric import Schwarzschild, Kerr, KerrNewman
from einsteinpy.geodesic import Geodesic, Timelike

from einsteinpy import constant

_c = constant.c.value
_G = constant.G.value
_Cc = constant.coulombs_const.value


def test_str_repr():
    """
    Tests, if the ``__str__`` and ``__repr__`` messages match

    """
    M = 1e25 * u.kg
    sph = SphericalDifferential(
        t=0.0 * u.s,
        r=306.0 * u.m,
        theta=np.pi / 2 * u.rad,
        phi=-np.pi / 2 * u.rad,
        v_r=0.0 * u.m / u.s,
        v_th=0.01 * u.rad / u.s,
        v_p=10.0 * u.rad / u.s,
    )
    ms = Schwarzschild(coords=sph, M=M)

    end_lambda = 1.
    step_size = 0.4e-6

    geod = Timelike(metric=ms, coords=sph, end_lambda=end_lambda, step_size=step_size)

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
    metric = Schwarzschild(coords=sph, M=6e24 * u.kg)

    end_lambda = 0.002
    step_size = 5e-8

    return sph, metric, end_lambda, step_size


def test_Geodesics_has_trajectory(dummy_data):
    sph, metric, end_lambda, step_size = dummy_data
    geo = Timelike(
        metric=metric,
        coords=sph,
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
            4e24 * u.kg,
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
            5.972e24 * u.kg,
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
            5.972e24 * u.kg,
            0.001,
            5e-6,
        ),
    ],
)
def test_calculate_trajectory_schwarzschild(sph, M, end_lambda, step_size):
    ms = Schwarzschild(coords=sph, M=M)

    geod = Timelike(
        metric=ms,
        coords=sph,
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
    testarray = np.array(testarray)

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
    metric = Schwarzschild(coords=sph, M=1.989e30 * u.kg)

    end_lambda = 3.154e7

    geod = Timelike(
        metric=metric,
        coords=sph,
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

    metric = Schwarzschild(coords=x_sph, M=1.989e30 * u.kg)

    end_lambda = 3.154e7

    geod = Timelike(
        metric=metric,
        coords=x_sph,
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
            4e24 * u.kg,
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
            5.972e24 * u.kg,
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
    metric = Schwarzschild(coords=sph, M=M)

    geod = Timelike(
        metric=metric,
        coords=sph,
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
    M = 1e25 * u.kg
    sph = SphericalDifferential(
        t=0.0 * u.s,
        r=306.0 * u.m,
        theta=np.pi / 2 * u.rad,
        phi=-np.pi / 2 * u.rad,
        v_r=0.0 * u.m / u.s,
        v_th=0.01 * u.rad / u.s,
        v_p=10.0 * u.rad / u.s,
    )
    ms = Schwarzschild(coords=sph, M=M)

    end_lambda = 1.
    stepsize = 0.4e-6
    OdeMethodKwargs = {"stepsize": stepsize}

    geod = Timelike(
        metric=ms,
        coords=sph,
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
            4e24 * u.kg,
            2e-3 * u.one,
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
            5.972e24 * u.kg,
            2e-3 * u.one,
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
            5.972e24 * u.kg,
            0. * u.one,
            0.001,
            5e-6,
        ),
    ],
)
def test_calculate_trajectory_kerr(
    bl, M, a, end_lambda, step_size
):
    mk = Kerr(coords=bl, M=M, a=a)

    geod = Timelike(
        metric=mk,
        coords=bl,
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
    testarray = np.array(testarray)

    assert_allclose(testarray, _c ** 2, 1e-8)


# Picking a test randomly and testing with Geodesic
# instead of Timelike
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
            4e24 * u.kg,
            2e-3 * u.one,
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
            5.972e24 * u.kg,
            2e-3 * u.one,
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
            5.972e24 * u.kg,
            0. * u.one,
            0.001,
            5e-6,
        ),
    ],
)
def test_calculate_trajectory_kerr_Geodesic(
    bl, M, a, end_lambda, step_size
):
    mk = Kerr(coords=bl, M=M, a=a)

    geod = Geodesic(
        time_like=True,
        metric=mk,
        coords=bl,
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
    testarray = np.array(testarray)

    assert_allclose(testarray, _c ** 2, 1e-8)


def test_calculate_trajectory3_kerr():
    # Based on the revolution of earth around sun
    # Data from https://en.wikipedia.org/wiki/Earth%27s_orbit
    # Initialized with cartesian coordinates
    # Function returning cartesian coordinates
    M = 1.989e30 * u.kg
    a = 0.0 * u.one

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

    mk = Kerr(coords=x_bl, M=M, a=a)

    end_lambda = 3.154e7

    geod = Timelike(
        metric=mk,
        coords=x_bl,
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
            4e24 * u.kg,
            2e-3 * u.one,
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
            5.972e24 * u.kg,
            0. * u.one,
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
    mk = Kerr(coords=bl, M=M, a=a)

    geod = Timelike(
        metric=mk,
        coords=bl,
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
    M = 1e25 * u.kg
    a = 0. * u.one

    bl = BoyerLindquistDifferential(
        t=0.0 * u.s,
        r=306.0 * u.m,
        theta=np.pi / 2 * u.rad,
        phi=np.pi / 2 * u.rad,
        v_r=0.0 * u.m / u.s,
        v_th=0.01 * u.rad / u.s,
        v_p=10.0 * u.rad / u.s,
    )

    mk = Kerr(coords=bl, M=M, a=a)

    end_lambda = 1.
    stepsize = 0.4e-6
    OdeMethodKwargs = {"stepsize": stepsize}

    geod = Timelike(
        metric=mk,
        coords=bl,
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
    M = 1.989e30 * u.kg
    a = 0. * u.one
    Q = 0. * u.C
    q = 0. * u.C / u.kg
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

    mkn = KerrNewman(coords=x_bl, M=M, a=a, Q=Q, q=q)

    end_lambda = 3.154e7

    geod = Timelike(
        metric=mkn,
        coords=x_bl,
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
    M = 0.5 * 5.972e24 * u.kg
    a = 0. * u.one
    Q = 11604461683.91822052001953125 * u.C
    q = _G * M.value / _Cc * u.C / u.kg

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

    mkn = KerrNewman(coords=bl, M=M, a=a, Q=Q, q=q)

    geod = Timelike(
        metric=mkn,
        coords=bl,
        end_lambda=end_lambda,
        step_size=step_size,
        return_cartesian=False
    )

    ans = geod.trajectory

    assert_allclose(ans[0][1], ans[-1][1], 1e-2)


@pytest.fixture()
def test_input():
    a = 1e-6 * u.one
    Q = 100. * u.C
    q = 1. * u.C / u.kg
    end_lambda = 200.0
    step_size = 1.0

    return a, Q, q, end_lambda, step_size


def test_compare_calculate_trajectory_iterator_bl_kerrnewman(test_input):
    a, Q, q, end_lambda, step_size = test_input
    M = 0.5 * 5.972e24 * u.kg

    bl = BoyerLindquistDifferential(
        t=0.0 * u.s,
        r=1e6 * u.m,
        theta=0.6 * np.pi * u.rad,
        phi=np.pi / 8 * u.rad,
        v_r=1e4 * u.m / u.s,
        v_th=-0.01 * u.rad / u.s,
        v_p=0.0 * u.rad / u.s,
    )

    mkn = KerrNewman(coords=bl, M=M, a=a, Q=Q, q=q)

    OdeMethodKwargs = {"stepsize": step_size}
    return_cartesian = False

    geod = Timelike(
        metric=mkn,
        coords=bl,
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
    M = 2e24 * u.kg

    x_bl = CartesianDifferential(
        t=0.0 * u.s,
        x=1e6 * u.m,
        y=1e6 * u.m,
        z=20.5 * u.m,
        v_x=1e4 * u.m / u.s,
        v_y=1e4 * u.m / u.s,
        v_z=-30.0 * u.m / u.s
    ).bl_differential(M=M, a=a)

    mkn = KerrNewman(coords=x_bl, M=M, a=a, Q=Q, q=q)

    OdeMethodKwargs = {"stepsize": step_size}
    return_cartesian = True

    geod = Timelike(
        metric=mkn,
        coords=x_bl,
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
    M = 0.5 * 5.972e24 * u.kg
    a = 0. * u.one
    Q = 0. * u.C
    q = 0. * u.C / u.kg

    bl = BoyerLindquistDifferential(
        t=0.0 * u.s,
        r=306.0 * u.m,
        theta=np.pi / 2 * u.rad,
        phi=np.pi / 2 * u.rad,
        v_r=0.0 * u.m / u.s,
        v_th=-0.01 * u.rad / u.s,
        v_p=10.0 * u.rad / u.s,
    )

    mkn = KerrNewman(coords=bl, M=M, a=a, Q=Q, q=q)

    end_lambda = 200.
    step_size = 0.4e-6
    OdeMethodKwargs = {"stepsize": step_size}

    geod = Timelike(metric=mkn, coords=bl, end_lambda=end_lambda, step_size=step_size)

    with warnings.catch_warnings(record=True) as w:
        it = geod.calculate_trajectory_iterator(
            OdeMethodKwargs=OdeMethodKwargs,
        )
        for _, _ in zip(range(1000), it):
            pass
        assert len(w) >= 1


def test_calculate_state_raises_TypeError():
    """
    Tests, if ``_calculate_state`` raises TypeError, in case of \
    coordinate mismatch

    """
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
    )  # .spherical_differential()

    metric = Schwarzschild(coords=x_sph.spherical_differential(), M=1.989e30 * u.kg)

    end_lambda = 3.154e7

    with pytest.raises(TypeError):
        geod = Timelike(
            metric=metric,
            coords=x_sph,
            end_lambda=end_lambda,
            step_size=end_lambda / 2e3,
        )
