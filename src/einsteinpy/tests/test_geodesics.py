import astropy.units as u
import numpy as np
import pytest
from numpy.testing import assert_allclose

from einsteinpy.bodies import Body
from einsteinpy.coordinates import SphericalDifferential
from einsteinpy.geodesic import Geodesic


@pytest.fixture()
def dummy_data():
    obj = SphericalDifferential(
        130 * u.m,
        np.pi / 2 * u.rad,
        -np.pi / 8 * u.rad,
        0 * u.m / u.s,
        0 * u.rad / u.s,
        1900 * u.rad / u.s,
    )
    att = Body(name="attractor", mass=6e24 * u.kg, parent=None)
    b1 = Body(name="obj", differential=obj, parent=att)
    t = 0 * u.s
    start_lambda = 0.0
    end_lambda = 0.002
    step_size = 5e-8
    return b1, t, start_lambda, end_lambda, step_size


def test_Geodesics_conserves_the_attractor(dummy_data):
    body, t, _, end_lambda, stepsize = dummy_data
    geo = Geodesic(body, time=t, end_lambda=end_lambda, step_size=stepsize)
    assert geo.attractor == body.parent


def test_Geodesics_has_trajectory(dummy_data):
    body, t, _, end_lambda, stepsize = dummy_data
    geo = Geodesic(body, time=t, end_lambda=end_lambda, step_size=stepsize)
    assert isinstance(geo.trajectory, np.ndarray)


def test_Geodesics_str_returns_members(dummy_data):
    body, t, _, end_lambda, stepsize = dummy_data
    geo = Geodesic(body, time=t, end_lambda=end_lambda, step_size=stepsize)
    assert (
        geo.__str__()
        == "Geodesic object body name= (obj) , metric=(Schwarzschild) , parent name=(attractor) , parent mass=(6e+24 kg)"
    )


def test_Geodesics_repr_returns_members(dummy_data):
    body, t, _, end_lambda, stepsize = dummy_data
    geo = Geodesic(body, time=t, end_lambda=end_lambda, step_size=stepsize)
    assert (
        geo.__repr__()
        == "Geodesic object body name= (obj) , metric=(Schwarzschild) , parent name=(attractor) , parent mass=(6e+24 kg)"
    )
