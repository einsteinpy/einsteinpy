from unittest import mock

import astropy.units as u
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pytest

from einsteinpy.plotting import StaticGeodesicPlotter


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


def test_staticgeodesicplotter_has_axes(dummy_data):
    r, v, _, m, _, el, ss = dummy_data
    cl = StaticGeodesicPlotter(m)
    assert isinstance(cl.ax, mpl.axes.SubplotBase)
    assert cl.time.value == 0.0
    assert cl._attractor_present is False


@mock.patch("einsteinpy.plotting.geodesics_static.StaticGeodesicPlotter.plot_attractor")
def test_plot_calls_plot_attractor(mock_plot_attractor):
    r = [306 * u.m, np.pi / 2 * u.rad, np.pi / 2 * u.rad]
    v = [0 * u.m / u.s, 0 * u.rad / u.s, 951.0 * u.rad / u.s]
    m = 4e24 * u.kg
    el = 0.002
    ss = 0.5e-6
    cl = StaticGeodesicPlotter(m)
    cl.plot(r, v, el, ss)
    mock_plot_attractor.assert_called_with()
    assert cl._attractor_present


@mock.patch("einsteinpy.plotting.geodesics_static.plt.savefig")
def test_plot_save_saves_plot(mock_save):
    r = [306 * u.m, np.pi / 2 * u.rad, np.pi / 2 * u.rad]
    v = [0 * u.m / u.s, 0 * u.rad / u.s, 951.0 * u.rad / u.s]
    m = 4e24 * u.kg
    el = 0.002
    ss = 0.5e-6
    cl = StaticGeodesicPlotter(m)
    cl.plot(r, v, el, ss)
    name = "test_plot.png"
    cl.save(name)
    mock_save.assert_called_with(name)
