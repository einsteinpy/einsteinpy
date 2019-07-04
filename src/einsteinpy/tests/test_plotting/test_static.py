from unittest import mock

import astropy.units as u
import numpy as np
import pytest
from matplotlib.axes import Axes

from einsteinpy.bodies import Body
from einsteinpy.coordinates import SphericalDifferential
from einsteinpy.geodesic import Geodesic
from einsteinpy.plotting import StaticGeodesicPlotter


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
    geo = Geodesic(b1, time=t, end_lambda=end_lambda, step_size=step_size)
    return geo


def test_staticgeodesicplotter_has_axes(dummy_data):
    geodesic = dummy_data
    cl = StaticGeodesicPlotter()
    assert isinstance(cl.ax, Axes)
    assert cl.attractor_present is False


@mock.patch("einsteinpy.plotting.geodesics.static.plt.show")
def test_plot_calls_plt_show(mock_show, dummy_data):
    geodesic = dummy_data
    cl = StaticGeodesicPlotter()
    cl.plot(geodesic)
    cl.show()
    mock_show.assert_called_with()
    assert cl.attractor_present


@mock.patch("einsteinpy.plotting.geodesics.static.plt.savefig")
def test_plot_save_saves_plot(mock_save, dummy_data):
    geodesic = dummy_data
    cl = StaticGeodesicPlotter()
    cl.plot(geodesic)
    name = "test_plot.png"
    cl.save(name)
    mock_save.assert_called_with(name)


def test_animate_creates_ani(dummy_data):
    geodesic = dummy_data
    cl = StaticGeodesicPlotter()
    cl.animate(geodesic, interval=10)
    assert cl.attractor_present
    assert cl.ani


def test_plot_calls_draw_attractor_Manualscale(dummy_data):
    geodesic = dummy_data
    cl = StaticGeodesicPlotter(attractor_radius_scale=1500)
    cl.plot(geodesic)
    assert cl.attractor_present
    assert cl.attractor_radius_scale == 1500


def test_plot_calls_draw_attractor_AutoScale(dummy_data):
    geodesic = dummy_data
    cl = StaticGeodesicPlotter()
    cl.plot(geodesic)
    assert cl.attractor_present
