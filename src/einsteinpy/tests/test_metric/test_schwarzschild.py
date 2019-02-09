import numpy as np
import pytest
from astropy import units as u
from numpy.testing import assert_allclose

from einsteinpy import constant
from einsteinpy.metric import Schwarzschild

_c = constant.c.value


@pytest.mark.parametrize(
    "pos_vec, vel_vec, time, M, start_lambda, end_lambda, OdeMethodKwargs",
    [
        (
            [306 * u.m, np.pi / 2 * u.rad, np.pi / 2 * u.rad],
            [0 * u.m / u.s, 0 * u.rad / u.s, 951.0 * u.rad / u.s],
            0 * u.s,
            4e24 * u.kg,
            0.0,
            0.002,
            {"stepsize": 0.5e-6},
        ),
        (
            [1 * u.km, 0.15 * u.rad, np.pi / 2 * u.rad],
            [_c * u.m / u.s, 0.5e-5 * _c * u.rad / u.s, 1e-4 * _c * u.rad / u.s],
            0 * u.s,
            5.972e24 * u.kg,
            0.0,
            0.000001,
            {"stepsize": 0.5e-6},
        ),
        (
            [50 * u.km, np.pi / 2 * u.rad, np.pi / 2 * u.rad],
            [0.1 * _c * u.m / u.s, 2e-7 * _c * u.rad / u.s, 1e-5 * u.rad / u.s],
            0 * u.s,
            5.972e24 * u.g,
            0.0,
            0.001,
            {"stepsize": 5e-6},
        ),
        (
            [1 * u.km, 0.8 * np.pi * u.rad, np.pi / 2 * u.rad],
            [1.3 * _c * u.m / u.s, 1e-6 * _c * u.rad / u.s, 3e-5 * _c * u.rad / u.s],
            0 * u.s,
            5.972e24 * u.kg,
            0.0,
            0.00001,
            {"stepsize": 0.5e-6},
        ),
    ],
)
def test_calculate_trajectory(
    pos_vec, vel_vec, time, M, start_lambda, end_lambda, OdeMethodKwargs
):
    cl = Schwarzschild.from_spherical(pos_vec, vel_vec, time, M)
    ans = cl.calculate_trajectory(
        start_lambda=start_lambda,
        end_lambda=end_lambda,
        OdeMethodKwargs=OdeMethodKwargs,
    )
    ans = ans[1]
    testarray = ans[:, 4] ** 2 - (
        (ans[:, 5] ** 2 + ans[:, 6] ** 2 + ans[:, 7] ** 2) / (_c ** 2)
    )
    comparearray = np.ones(shape=ans[:, 4].shape, dtype=float)
    assert_allclose(testarray, comparearray, 1e-3)


def test_calculate_trajectory2():
    # based on the revolution of earth around sun
    # data from https://en.wikipedia.org/wiki/Earth%27s_orbit
    M = 1.989e30 * u.kg
    distance_at_perihelion = 147.10e6 * u.km
    speed_at_perihelion = 30.29 * u.km / u.s
    angular_vel = (speed_at_perihelion / distance_at_perihelion) * u.rad
    pos_vec = [distance_at_perihelion, np.pi / 2 * u.rad, 0 * u.rad]
    vel_vec = [0 * u.km / u.s, 0 * u.rad / u.s, angular_vel]
    end_lambda = ((1 * u.year).to(u.s)).value
    cl = Schwarzschild.from_spherical(pos_vec, vel_vec, 0 * u.s, M)
    ans = cl.calculate_trajectory(
        start_lambda=0.0,
        end_lambda=end_lambda,
        OdeMethodKwargs={"stepsize": end_lambda / 5e3},
    )[1]
    # velocity should be 29.29 km/s at apehelion(where r is max)
    i = np.argmax(ans[:, 1])  # index whre radial distance is max
    v_apehelion = (((ans[i][1] * ans[i][7]) * (u.m / u.s)).to(u.km / u.s)).value
    assert_allclose(v_apehelion, 29.29, rtol=0.01)


def test_calculate_trajectory3():
    # same test as with test_calculate_trajectory2(),
    # but initialialized with cartesian coordinates
    M = 1.989e30 * u.kg
    distance_at_perihelion = 147.10e6 * u.km
    speed_at_perihelion = 30.29 * u.km / u.s
    pos_vec = [
        distance_at_perihelion / np.sqrt(2),
        distance_at_perihelion / np.sqrt(2),
        0 * u.km,
    ]
    vel_vec = [
        -1 * speed_at_perihelion / np.sqrt(2),
        speed_at_perihelion / np.sqrt(2),
        0 * u.km / u.h,
    ]
    end_lambda = ((1 * u.year).to(u.s)).value
    cl = Schwarzschild.from_cartesian(pos_vec, vel_vec, 0 * u.min, M)
    ans = cl.calculate_trajectory(
        start_lambda=0.0,
        end_lambda=end_lambda,
        OdeMethodKwargs={"stepsize": end_lambda / 5e3},
    )[1]
    # velocity should be 29.29 km/s at apehelion(where r is max)
    i = np.argmax(ans[:, 1])  # index whre radial distance is max
    v_apehelion = (((ans[i][1] * ans[i][7]) * (u.m / u.s)).to(u.km / u.s)).value
    assert_allclose(v_apehelion, 29.29, rtol=0.01)
