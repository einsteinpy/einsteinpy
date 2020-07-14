import astropy.units as u
import numpy as np

from numpy.testing import assert_allclose

from einsteinpy.units import primitive


def test_primitive():
    """
    Tests, if quantities are safely being converted into ``numpy.float64`` values

    """
    M = 10 * u.kg
    a = 0.29 * u.one
    Q = 1e3 * u.C
    t = 1e3 * u.s
    r = 1e10 * u.m
    theta = np.pi / 8 * u.rad
    temp = 1e5

    values = [10., 0.29, 1e3, 1e3, 1e10, np.pi / 8, 1e5]
    values_calc = primitive(M, a, Q, t, r, theta, temp)

    assert_allclose(values, values_calc)
