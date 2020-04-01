import warnings

import numpy as np
import pytest
from astropy import units as u
from numpy.testing import assert_allclose

from einsteinpy import constant
from einsteinpy.coordinates import BoyerLindquistDifferential, CartesianDifferential
from einsteinpy.metric import Kerr
from einsteinpy.utils import kerr_utils, schwarzschild_radius_dimensionless

_c = constant.c.value


@pytest.mark.parametrize(
    "coords, time, M, start_lambda, end_lambda, OdeMethodKwargs",
    [
        (
            BoyerLindquistDifferential(
                306 * u.m,
                np.pi / 2.05 * u.rad,
                np.pi / 2 * u.rad,
                0 * u.m / u.s,
                0 * u.rad / u.s,
                951.0 * u.rad / u.s,
                2e-3 * u.m,
            ),
            0 * u.s,
            4e24 * u.kg,
            0.0,
            0.001,
            {"stepsize": 0.5e-6},
        ),
        (
            BoyerLindquistDifferential(
                1 * u.km,
                0.15 * u.rad,
                np.pi / 2 * u.rad,
                0.1 * _c * u.m / u.s,
                0.5e-5 * _c * u.rad / u.s,
                0.5e-4 * _c * u.rad / u.s,
                2e-3 * u.m,
            ),
            0 * u.s,
            5.972e24 * u.kg,
            0.0,
            0.0001,
            {"stepsize": 0.5e-6},
        ),
        (
            BoyerLindquistDifferential(
                50 * u.km,
                np.pi / 2 * u.rad,
                np.pi / 2 * u.rad,
                0.1 * _c * u.m / u.s,
                2e-7 * _c * u.rad / u.s,
                1e-5 * u.rad / u.s,
                0.0 * u.km,
            ),
            0 * u.s,
            5.972e24 * u.g,
            0.0,
            0.001,
            {"stepsize": 5e-6},
        ),
    ],
)
def test_calculate_trajectory(
    coords, time, M, start_lambda, end_lambda, OdeMethodKwargs
):
    _scr = schwarzschild_radius_dimensionless(M)
    obj = Kerr.from_coords(coords, M, time)
    ans = obj.calculate_trajectory(
        start_lambda=start_lambda,
        end_lambda=end_lambda,
        OdeMethodKwargs=OdeMethodKwargs,
    )
    ans = ans[1]
    testarray = list()
    for i in ans:
        g = kerr_utils.metric(i[1], i[2], M.value, coords.a.to(u.m).value)
        testarray.append(
            g[0][0] * (i[4] ** 2)
            + g[1][1] * (i[5] ** 2)
            + g[2][2] * (i[6] ** 2)
            + g[3][3] * (i[7] ** 2)
            + 2 * g[0][3] * i[4] * i[7]
        )
    testarray = np.array(testarray, dtype=float)
    comparearray = np.ones(shape=ans[:, 4].shape, dtype=float)
    assert_allclose(testarray, comparearray, 1e-4)


def test_calculate_trajectory3():
    # Based on the revolution of earth around sun
    # Data from https://en.wikipedia.org/wiki/Earth%27s_orbit
    # Initialialized with cartesian coordinates
    # Function returning cartesian coordinates
    M = 1.989e30 * u.kg
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
    end_lambda = ((1 * u.year).to(u.s)).value
    a = 0 * u.km
    cl = Kerr.from_coords(cart_obj, M, a)
    ans = cl.calculate_trajectory(
        start_lambda=0.0,
        end_lambda=end_lambda,
        return_cartesian=True,
        OdeMethodKwargs={"stepsize": end_lambda / 2e3},
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
    "coords, time, M, start_lambda, end_lambda, OdeMethodKwargs, return_cartesian",
    [
        (
            BoyerLindquistDifferential(
                306 * u.m,
                np.pi / 2 * u.rad,
                np.pi / 2 * u.rad,
                0 * u.m / u.s,
                0.1 * u.rad / u.s,
                951.0 * u.rad / u.s,
                2e-3 * u.m,
            ),
            0 * u.s,
            4e24 * u.kg,
            0.0,
            0.0003,
            {"stepsize": 0.3e-6},
            True,
        ),
        (
            BoyerLindquistDifferential(
                1 * u.km,
                0.15 * u.rad,
                np.pi / 2 * u.rad,
                0.2 * _c * u.m / u.s,
                0.5e-5 * _c * u.rad / u.s,
                1e-4 * _c * u.rad / u.s,
                0.0 * u.km,
            ),
            0 * u.s,
            5.972e24 * u.kg,
            0.0,
            0.0004,
            {"stepsize": 0.5e-6},
            False,
        ),
    ],
)
def test_calculate_trajectory_iterator(
    coords, time, M, start_lambda, end_lambda, OdeMethodKwargs, return_cartesian
):
    cl1 = Kerr.from_coords(coords, M, time)
    arr1 = cl1.calculate_trajectory(
        start_lambda=start_lambda,
        end_lambda=end_lambda,
        OdeMethodKwargs=OdeMethodKwargs,
        return_cartesian=return_cartesian,
    )[1]
    cl2 = Kerr.from_coords(coords, M, time)
    it = cl2.calculate_trajectory_iterator(
        start_lambda=start_lambda,
        OdeMethodKwargs=OdeMethodKwargs,
        return_cartesian=return_cartesian,
    )
    arr2_list = list()
    for _, val in zip(range(100), it):
        arr2_list.append(val[1])
    arr2 = np.array(arr2_list)
    assert_allclose(arr1[:100, :], arr2, rtol=1e-10)


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
    OdeMethodKwargs = {"stepsize": 0.4e-6}
    cl = Kerr.from_coords(bl_obj, M)
    with warnings.catch_warnings(record=True) as w:
        it = cl.calculate_trajectory_iterator(
            start_lambda=start_lambda,
            OdeMethodKwargs=OdeMethodKwargs,
            stop_on_singularity=True,
        )
        for _, _ in zip(range(1000), it):
            pass
        assert len(w) >= 1
