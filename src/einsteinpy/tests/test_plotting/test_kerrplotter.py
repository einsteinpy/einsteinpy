from unittest import mock

import astropy.units as u
import numpy as np
import pytest

from einsteinpy.plotting import KerrPlotter


@mock.patch("einsteinpy.plotting.kerr_plot.KerrPlotter._calc_event_horizon")
def test_plot_event_horizon_calls_calc_event_horizon(mock_calc_event_horizon):
    m = 4e30
    start = 0
    end = np.pi
    steps = 1000
    color = "blue"
    coord = "Spherical"
    obj = KerrPlotter(m)
    obj.plot_event_horizon(start, end, steps, color, coord)
    mock_calc_event_horizon.assert_called_with(start, end, steps, coord)


@mock.patch("einsteinpy.plotting.kerr_plot.KerrPlotter._calc_ergosphere")
def test_plot_ergosphere_calls_calc_ergosphere(mock_calc_ergosphere):
    m = 4e30
    start = 0
    end = np.pi
    steps = 1000
    color = "red"
    coord = "Spherical"
    obj = KerrPlotter(m)
    obj.plot_ergosphere(start, end, steps, color, coord)
    mock_calc_ergosphere.assert_called_with(start, end, steps, coord)


@mock.patch("einsteinpy.plotting.kerr_plot.KerrPlotter._calc_ergosphere")
@mock.patch("einsteinpy.plotting.kerr_plot.KerrPlotter._calc_event_horizon")
def test_plot_calls_calc_functions(mock_calc_event_horizon, mock_calc_ergosphere):
    m = 4e30
    start = 0
    end = np.pi
    steps = 1000
    e_color = "red"
    h_color = "blue"
    coord = "Spherical"
    obj = KerrPlotter(m)
    obj.plot(start, end, steps, e_color, h_color, coord)
    mock_calc_ergosphere.assert_called_with(start, end, steps, coord)
    mock_calc_event_horizon.assert_called_with(start, end, steps, coord)


@mock.patch("einsteinpy.plotting.kerr_plot.plt.show")
def test_plot_show_shows_plot(mock_show):
    m = 4e30
    start = 0
    end = np.pi
    steps = 1000
    e_color = "red"
    h_color = "blue"
    obj = KerrPlotter(m)
    obj.plot(start, end, steps, e_color, h_color)
    obj.show()
    mock_show.assert_called_with()
