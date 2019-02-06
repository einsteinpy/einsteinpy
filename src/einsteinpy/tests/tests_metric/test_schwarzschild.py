import pytest
from astropy import units as u
from astropy.tests.helper import assert_quantity_allclose

from einsteinpy.metric import schwarzschild


def test_constructor_values():
    # pos_vec = # in astropy units
    # vel_vec =
    # time =
    # M =
    # metric_obj = Schwarzschild(pos_vec, vel_vec, time, M)
    # assert metric_obj.pos_vec == pos_vec
    # assert metric_obj.vel_vec == vel_vec
    # assert metric_obj.time == time
    # assert metric_obj.M == M
    # ### Uncomment when completely ready
    pass


def test_christoffel_symbols_values():
    pass


def test_f_values():
    pass


def test_f_vec_values():
    pass


def test_calculate_trajectory_four_velocity_constant():
    pass


def test_calculate_trajectory_vec_values():
    pass
