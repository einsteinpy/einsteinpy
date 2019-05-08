import warnings

import numpy as np
import pytest
from astropy import units as u
from einsteinpy import constant
from einsteinpy.metric import KerrNewman
from einsteinpy.utils import kerrnewman_utils, schwarzschild_radius
from numpy.testing import assert_allclose

_c = constant.c.value
_G = constant.G.value
_cc = constant.coulombs_const.value


def test_calculate_trajectory0():
    # Based on the revolution of earth around sun
    # Data from https://en.wikipedia.org/wiki/Earth%27s_orbit
    # Initialialized with cartesian coordinates
    # Function returning cartesian coordinates
    M = 1.989e30 * u.kg
    q = 0 * u.C / u.kg
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
    cl = KerrNewman.from_cartesian(pos_vec, vel_vec, q, 0 * u.min, M, 0.0, 0 * u.C)
    ans = cl.calculate_trajectory(
        start_lambda=0.0,
        end_lambda=end_lambda,
        return_cartesian=True,
        OdeMethodKwargs={"stepsize": end_lambda / 1.5e3},
    )[1]
    # velocity should be 29.29 km/s at apehelion(where r is max)
    R = np.sqrt(ans[:, 1] ** 2 + ans[:, 2] ** 2 + ans[:, 3] ** 2)
    i = np.argmax(R)  # index whre radial distance is max
    v_apehelion = (
        (np.sqrt(ans[i, 5] ** 2 + ans[i, 6] ** 2 + ans[i, 7] ** 2) * (u.m / u.s)).to(
            u.km / u.s
        )
    ).value
    assert_allclose(v_apehelion, 29.29, rtol=0.01)


@pytest.mark.parametrize(
    "M, r, end_lambda, stepsize", [(0.5 * 5.972e24 * u.kg, 1000.0 * u.km, 1000.0, 0.5)]
)
def test_calculate_trajectory1(M, r, end_lambda, stepsize):
    # the test particle should not move as gravitational & electromagnetic forces are balanced
    tmp = _G * M.value / _cc
    q = tmp * u.C / u.kg
    print(tmp, 5.900455 * tmp ** 3)
    Q = 11604461683.91822052001953125 * u.C
    pos_vec = [r, 0.5 * np.pi * u.rad, 0 * u.rad]
    vel_vec = [0 * u.m / u.s, 0 * u.rad / u.s, 0.0 * u.rad / u.s]
    cl = KerrNewman.from_BL(pos_vec, vel_vec, q, 0 * u.s, M, 0.0, Q)
    ans = cl.calculate_trajectory(
        end_lambda=end_lambda, OdeMethodKwargs={"stepsize": stepsize}
    )
    assert_allclose(ans[1][0][1], ans[1][-1][1], 1e-2)


@pytest.mark.parametrize(
    "pos_vec, vel_vec, q, M, a, Q, el, ss",
    [
        (
            [1000.0 * u.km, 0.6 * np.pi * u.rad, np.pi / 8 * u.rad],
            [10000 * u.m / u.s, -0.01 * u.rad / u.s, 0.0 * u.rad / u.s],
            1 * u.C / u.g,
            0.5 * 5.972e24 * u.kg,
            1e-6,
            100 * u.C,
            200.0,
            1.0,
        )
    ],
)
def test_compare_calculate_trajectory_iterator_bl(pos_vec, vel_vec, q, M, a, Q, el, ss):
    cl1 = KerrNewman.from_BL(pos_vec, vel_vec, q, 0 * u.s, M, a, Q)
    cl2 = KerrNewman.from_BL(pos_vec, vel_vec, q, 0 * u.s, M, a, Q)
    ans1 = cl1.calculate_trajectory(end_lambda=el, OdeMethodKwargs={"stepsize": ss})[1]
    it = cl2.calculate_trajectory_iterator(OdeMethodKwargs={"stepsize": ss})
    ans2 = list()
    for _, val in zip(range(20), it):
        ans2.append(val[1])
    ans2 = np.array(ans2)
    print(ans1)
    assert_allclose(ans1[:20], ans2)


@pytest.mark.parametrize(
    "pos_vec, vel_vec, q, M, a, Q, el, ss",
    [
        (
            [1000000 * u.m, 1000000 * u.m, 20.5 * u.m],
            [10000 * u.m / u.s, 10000 * u.m / u.s, -30 * u.m / u.s],
            1 * u.C / u.g,
            2e24 * u.kg,
            1e-6,
            100 * u.C,
            200.0,
            1.0,
        )
    ],
)
def test_compare_calculate_trajectory_iterator_cartesians(
    pos_vec, vel_vec, q, M, a, Q, el, ss
):
    cl1 = KerrNewman.from_cartesian(pos_vec, vel_vec, q, 0 * u.s, M, a, Q)
    cl2 = KerrNewman.from_cartesian(pos_vec, vel_vec, q, 0 * u.s, M, a, Q)
    ans1 = cl1.calculate_trajectory(
        end_lambda=el, OdeMethodKwargs={"stepsize": ss}, return_cartesian=True
    )[1]
    it = cl2.calculate_trajectory_iterator(
        OdeMethodKwargs={"stepsize": ss}, return_cartesian=True
    )
    ans2 = list()
    for _, val in zip(range(20), it):
        ans2.append(val[1])
    ans2 = np.array(ans2)
    print(ans1)
    assert_allclose(ans1[:20], ans2)
