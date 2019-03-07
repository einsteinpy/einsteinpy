from unittest import mock

import astropy.units as u
import numpy as np
import pytest

from einsteinpy.plotting import ScatterGeodesicPlotter


@pytest.fixture()
def dummy_data():
    r = [306 * u.m, np.pi / 2 * u.rad, np.pi / 2 * u.rad]
    v = [0 * u.m / u.s, 0 * u.rad / u.s, 951.0 * u.rad / u.s]
    t = 0 * u.s
    m = 4e24 * u.kg
    start_lambda = 0.0
    end_lambda = 0.002
    step_size = 0.5e-6
    return r, v, t, m, start_lambda, end_lambda, step_size


def test_plot_attractor_is_called_only_once(dummy_data):
    r, v, _, m, _, el, ss = dummy_data
    cl = ScatterGeodesicPlotter(m)
    assert cl._attractor_present == False
    cl.plot(r, v, el, ss)
    assert cl._attractor_present == True


@mock.patch(
    "einsteinpy.plotting.geodesics_static.ScatterGeodesicPlotter._plot_attractor"
)
def test_plot_calls_plot_attractor(mock_plot_attractor):
    r = [306 * u.m, np.pi / 2 * u.rad, np.pi / 2 * u.rad]
    v = [0 * u.m / u.s, 0 * u.rad / u.s, 951.0 * u.rad / u.s]
    m = 4e24 * u.kg
    el = 0.002
    ss = 0.5e-6
    cl = ScatterGeodesicPlotter(m)
    cl.plot(r, v, el, ss)
    mock_plot_attractor.assert_called_with()


@mock.patch("einsteinpy.plotting.geodesics_static.plt.show")
def test_plot_show_shows_plot(mock_show):
    r = [306 * u.m, np.pi / 2 * u.rad, np.pi / 2 * u.rad]
    v = [0 * u.m / u.s, 0 * u.rad / u.s, 951.0 * u.rad / u.s]
    m = 4e24 * u.kg
    el = 0.002
    ss = 0.5e-6
    cl = ScatterGeodesicPlotter(m)
    cl.plot(r, v, el, ss)
    cl.show()
    mock_show.assert_called_with()


@mock.patch("einsteinpy.plotting.geodesics_static.plt.savefig")
def test_plot_save_saves_plot(mock_save):
    r = [306 * u.m, np.pi / 2 * u.rad, np.pi / 2 * u.rad]
    v = [0 * u.m / u.s, 0 * u.rad / u.s, 951.0 * u.rad / u.s]
    m = 4e24 * u.kg
    el = 0.002
    ss = 0.5e-6
    cl = ScatterGeodesicPlotter(m)
    cl.plot(r, v, el, ss)
    name = "test_plot.png"
    cl.save(name)
    mock_save.assert_called_with(name)
