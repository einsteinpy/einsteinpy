import numpy as np
import pytest
from astropy import units as u
from numpy.testing import assert_allclose

from einsteinpy import constant
from einsteinpy.metric import RobertsonWalker

_c = constant.c.value

def test_calculate_trajectory(pos_vec, vel_vec, time, M, era, tuning_param, k, start_lambda, end_lambda, OdeMethodKwargs):
    cl = RobertsonWalker.from_spherical(pos_vec, vel_vec, time, M, era, tuning_param, k)
    ans = cl.calculate_trajectory(start_lambda=start_lambda, end_lambda=end_lambda, OdeMethodKwargs=OdeMethodKwargs)
    ans=ans[1]
    testarray = ans[:, 4] ** 2 - (
        (ans[:, 5] ** 2 + ans[:, 6] ** 2 + ans[:, 7] ** 2) / (_c ** 2)
    )
    print(testarray)

test_calculate_trajectory([306 * u.m, np.pi / 2 * u.rad, np.pi / 2 * u.rad], [0 * u.m / u.s, 0 * u.rad / u.s, 951.0 * u.rad / u.s], 10000*u.s, 1*u.g, 'md', 
    1., 0, 0.0, 0.1, {'stepsize':0.0001})