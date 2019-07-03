from unittest import mock

import astropy.units as u
import numpy as np
import pytest

from einsteinpy.bodies import Body
from einsteinpy.coordinates import SphericalDifferential
from einsteinpy.geodesic import Geodesic
from einsteinpy.plotting.senile import ScatterGeodesicPlotter


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


def test_plot_attractor_is_called_only_once(dummy_data):
    geodesics = dummy_data
    cl = ScatterGeodesicPlotter()
    assert not cl._attractor_present
    cl.plot(geodesics)
    assert cl._attractor_present


@mock.patch(
    "einsteinpy.plotting.senile.geodesics_scatter.ScatterGeodesicPlotter._plot_attractor"
)
def test_plot_calls_plot_attractor(mock_plot_attractor, dummy_data):
    geodesics = dummy_data
    cl = ScatterGeodesicPlotter()
    cl.plot(geodesics)
    mock_plot_attractor.assert_called_with()


@mock.patch("einsteinpy.plotting.senile.geodesics_scatter.plt.figure")
def test_animate_calls_figure(mock_figure, dummy_data):
    geodesics = dummy_data
    cl = ScatterGeodesicPlotter()
    cl.animate(geodesics)
    mock_figure.assert_called_with()


@mock.patch("einsteinpy.plotting.senile.geodesics_scatter.plt.show")
def test_plot_show_shows_plot(mock_show, dummy_data):
    geodesics = dummy_data
    cl = ScatterGeodesicPlotter()
    cl.plot(geodesics)
    cl.show()
    mock_show.assert_called_with()


@mock.patch("einsteinpy.plotting.senile.geodesics_scatter.plt.savefig")
def test_plot_save_saves_plot(mock_save, dummy_data):
    geodesics = dummy_data
    cl = ScatterGeodesicPlotter()
    cl.plot(geodesics)
    name = "test_plot.png"
    cl.save(name)
    mock_save.assert_called_with(name)
