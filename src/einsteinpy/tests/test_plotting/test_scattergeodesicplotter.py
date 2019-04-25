from unittest import mock

import astropy.units as u
import numpy as np
import pytest

from einsteinpy.coordinates import SphericalDifferential
from einsteinpy.plotting import ScatterGeodesicPlotter


@pytest.fixture()
def dummy_data():
    sph_obj = SphericalDifferential(
        306 * u.m,
        np.pi / 2 * u.rad,
        np.pi / 2 * u.rad,
        0 * u.m / u.s,
        0 * u.rad / u.s,
        951.0 * u.rad / u.s,
    )
    t = 0 * u.s
    m = 4e24 * u.kg
    start_lambda = 0.0
    end_lambda = 0.002
    step_size = 0.5e-6
    return sph_obj, t, m, start_lambda, end_lambda, step_size


def test_plot_attractor_is_called_only_once(dummy_data):
    sph_obj, _, m, _, el, ss = dummy_data
    cl = ScatterGeodesicPlotter(m)
    assert not cl._attractor_present
    cl.plot(sph_obj, el, ss)
    assert cl._attractor_present


@mock.patch(
    "einsteinpy.plotting.geodesics_scatter.ScatterGeodesicPlotter._plot_attractor"
)
def test_plot_calls_plot_attractor(mock_plot_attractor, dummy_data):
    sph_obj, _, m, _, el, ss = dummy_data
    cl = ScatterGeodesicPlotter(m)
    cl.plot(sph_obj, el, ss)
    mock_plot_attractor.assert_called_with()


@mock.patch("einsteinpy.plotting.geodesics_scatter.plt.figure")
def test_animate_calls_figure(mock_figure, dummy_data):
    sph_obj, _, m, _, el, ss = dummy_data
    cl = ScatterGeodesicPlotter(m)
    cl.animate(sph_obj, el, ss)
    mock_figure.assert_called_with()


@mock.patch("einsteinpy.plotting.geodesics_scatter.plt.show")
def test_plot_show_shows_plot(mock_show, dummy_data):
    sph_obj, _, m, _, el, ss = dummy_data
    cl = ScatterGeodesicPlotter(m)
    cl.plot(sph_obj, el, ss)
    cl.show()
    mock_show.assert_called_with()


@mock.patch("einsteinpy.plotting.geodesics_scatter.plt.savefig")
def test_plot_save_saves_plot(mock_save, dummy_data):
    sph_obj, _, m, _, el, ss = dummy_data
    cl = ScatterGeodesicPlotter(m)
    cl.plot(sph_obj, el, ss)
    name = "test_plot.png"
    cl.save(name)
    mock_save.assert_called_with(name)
