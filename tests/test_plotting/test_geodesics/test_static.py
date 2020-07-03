from unittest import mock

import astropy.units as u
import numpy as np
import pytest
from matplotlib.axes import Axes

from einsteinpy.coordinates import SphericalDifferential
from einsteinpy.metric import Schwarzschild
from einsteinpy.geodesic import TimelikeGeodesic
from einsteinpy.plotting import StaticGeodesicPlotter


@pytest.fixture()
def dummy_data():
    M = 6e24 * u.kg

    sph = SphericalDifferential(
        t=0.0 * u.s,
        r=130.0 * u.m,
        theta=np.pi / 2 * u.rad,
        phi=-np.pi / 8 * u.rad,
        v_r=0.0 * u.m / u.s,
        v_th=0.0 * u.rad / u.s,
        v_p=1900.0 * u.rad / u.s,
    )

    ms = Schwarzschild(coords=sph, M=M)

    end_lambda = 0.002
    step_size = 5e-8

    geod = TimelikeGeodesic(metric=ms, coords=sph, end_lambda=end_lambda, step_size=step_size)

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
