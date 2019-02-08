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
    print(max(abs(comparearray - testarray)))
    assert_allclose(testarray, comparearray, 1e-3)


# def test_constructor_values():
#     pass


# def test_christoffel_symbols_values():
#     pass


# def test_f_values():
#     pass


# def test_f_vec_values():
#     pass


# def test_calculate_trajectory_four_velocity_constant():
#     pass


# def test_calculate_trajectory_vec_values():
#     pass
