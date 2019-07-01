import warnings

import numpy as np
import pytest
from astropy import units as u
from numpy.testing import assert_allclose

from einsteinpy import constant
from einsteinpy.coordinates import BoyerLindquistDifferential, CartesianDifferential
from einsteinpy.metric import KerrNewman
from einsteinpy.utils import kerrnewman_utils, schwarzschild_radius

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
    Q = 0 * u.C
    distance_at_perihelion = 147.10e6 * u.km
    speed_at_perihelion = 30.29 * u.km / u.s
    cart_obj = CartesianDifferential(
        distance_at_perihelion / np.sqrt(2),
        distance_at_perihelion / np.sqrt(2),
        0 * u.km,
        -1 * speed_at_perihelion / np.sqrt(2),
        speed_at_perihelion / np.sqrt(2),
        0 * u.km / u.h,
    )
    a = 0 * u.m
    end_lambda = ((1 * u.year).to(u.s)).value
    cl = KerrNewman.from_coords(coords=cart_obj, M=M, q=q, a=a, Q=Q)
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


def test_calculate_trajectory1():
    # the test particle should not move as gravitational & electromagnetic forces are balanced
    M = 0.5 * 5.972e24 * u.kg
    r = 1000.0 * u.km
    end_lambda = 1000.0
    stepsize = 0.5
    tmp = _G * M.value / _cc
    q = tmp * u.C / u.kg
    a = 0 * u.km
    Q = 11604461683.91822052001953125 * u.C
    bl_obj = BoyerLindquistDifferential(
        r,
        0.5 * np.pi * u.rad,
        0 * u.rad,
        0 * u.m / u.s,
        0 * u.rad / u.s,
        0.0 * u.rad / u.s,
        a,
    )
    cl = KerrNewman.from_coords(coords=bl_obj, q=q, M=M, Q=Q)
    ans = cl.calculate_trajectory(
        end_lambda=end_lambda, OdeMethodKwargs={"stepsize": stepsize}
    )
    assert_allclose(ans[1][0][1], ans[1][-1][1], 1e-2)


@pytest.fixture()
def test_input():
    q = 1 * u.C / u.g
    a = 1e-6 * u.m
    Q = 100 * u.C
    el = 200.0
    ss = 1.0
    return q, a, Q, el, ss


def test_compare_calculate_trajectory_iterator_bl(test_input):
    q, a, Q, el, ss = test_input
    bl_obj = BoyerLindquistDifferential(
        1000.0 * u.km,
        0.6 * np.pi * u.rad,
        np.pi / 8 * u.rad,
        10000 * u.m / u.s,
        -0.01 * u.rad / u.s,
        0.0 * u.rad / u.s,
        a,
    )
    M = 0.5 * 5.972e24 * u.kg
    cl1 = KerrNewman.from_coords(coords=bl_obj, q=q, M=M, Q=Q)
    cl2 = KerrNewman.from_coords(coords=bl_obj, q=q, M=M, Q=Q)
    ans1 = cl1.calculate_trajectory(end_lambda=el, OdeMethodKwargs={"stepsize": ss})[1]
    it = cl2.calculate_trajectory_iterator(OdeMethodKwargs={"stepsize": ss})
    ans2 = list()
    for _, val in zip(range(20), it):
        ans2.append(val[1])
    ans2 = np.array(ans2)
    assert_allclose(ans1[:20], ans2)


def test_compare_calculate_trajectory_iterator_cartesians(test_input):
    cart_obj = CartesianDifferential(
        1000000 * u.m,
        1000000 * u.m,
        20.5 * u.m,
        10000 * u.m / u.s,
        10000 * u.m / u.s,
        -30 * u.m / u.s,
    )
    M = 2e24 * u.kg
    q, a, Q, el, ss = test_input
    cl1 = KerrNewman.from_coords(coords=cart_obj, q=q, M=M, a=a, Q=Q)
    cl2 = KerrNewman.from_coords(coords=cart_obj, q=q, M=M, a=a, Q=Q)
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
    assert_allclose(ans1[:20], ans2)


def test_calculate_trajectory_iterator_RuntimeWarning():
    bl_obj = BoyerLindquistDifferential(
        306 * u.m,
        np.pi / 2 * u.rad,
        np.pi / 2 * u.rad,
        0 * u.m / u.s,
        0.01 * u.rad / u.s,
        10 * u.rad / u.s,
        0 * u.m,
    )
    M = 1e25 * u.kg
    start_lambda = 0.0
    q, Q = 0 * u.C / u.kg, 0 * u.C
    OdeMethodKwargs = {"stepsize": 0.4e-6}
    cl = KerrNewman.from_coords(coords=bl_obj, q=q, M=M, Q=Q)
    with warnings.catch_warnings(record=True) as w:
        it = cl.calculate_trajectory_iterator(
            start_lambda=start_lambda,
            OdeMethodKwargs=OdeMethodKwargs,
            stop_on_singularity=True,
        )
        for _, _ in zip(range(1000), it):
            pass
        assert len(w) >= 1
