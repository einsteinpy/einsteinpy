from unittest import mock

import numpy as np
import pytest
from matplotlib.axes import Axes

from einsteinpy.metric import Schwarzschild
from einsteinpy.coordinates.utils import four_position, stacked_vec
from einsteinpy.geodesic import Geodesic
from einsteinpy.plotting import StaticGeodesicPlotter


@pytest.fixture()
def dummy_data():
    M = 6e24
    t = 0.
    x_vec = np.array([130.0, np.pi / 2, -np.pi / 8])
    v_vec = np.array([0.0, 0.0, 1900.0])

    ms_cov = Schwarzschild(M=M)
    x_4vec = four_position(t, x_vec)
    ms_cov_mat = ms_cov.metric_covariant(x_4vec)
    init_vec = stacked_vec(ms_cov_mat, t, x_vec, v_vec, time_like=True)

    end_lambda = 0.002
    step_size = 5e-8

    geod = Geodesic(metric=ms_cov, init_vec=init_vec, end_lambda=end_lambda, step_size=step_size)

    return geod


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


@mock.patch("einsteinpy.plotting.geodesics.static.plt.show")
def test_plot_in_3d(mock_show, dummy_data):
    geodesic = dummy_data
    cl = StaticGeodesicPlotter(use_3d=True)
    cl.plot(geodesic)
    cl.show()
    mock_show.assert_called_with()
    assert cl.attractor_present
    assert cl.ax.name == "3d"


def test_set_scaling(dummy_data):
    geodesic = dummy_data
    cl = StaticGeodesicPlotter(use_3d=True)
    cl.plot(geodesic)
    cl._set_scaling(10, 1e-6, 10, 1e-5)
    cl._set_scaling(1e-6, 1e-6, 1e-6, 1e-5)
    zmin, zmax = cl.ax.get_zlim()
    ymin, ymax = cl.ax.get_ylim()
    xmin, xmax = cl.ax.get_xlim()
    assert xmin != -1e-5 and xmax != 1e-5
    assert ymin == -1e-5 and ymax == 1e-5
    assert zmin == -1e-5 and zmax == 1e-5


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
