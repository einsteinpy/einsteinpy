from unittest import mock

import numpy as np
import pytest
import warnings
from matplotlib.axes import Axes
from matplotlib import pyplot as plt

from einsteinpy.geodesic import Timelike
from einsteinpy.plotting import StaticGeodesicPlotter


@pytest.fixture()
def dummy_geod():
    return Timelike(
        metric="Kerr",
        metric_params=(0.9,),
        position=[2.15, np.pi / 2, 0.],
        momentum=[0., 0., 1.5],
        steps=50,
        delta=0.5,
        omega=0.01,  # Close orbit
        return_cartesian=True,
        suppress_warnings=True,
    )


def test_static_geod_plotter_has_axes(dummy_geod):
    geod = dummy_geod

    sgpl = StaticGeodesicPlotter()
    sgpl.plot(geod)

    assert isinstance(sgpl.ax, Axes)


@mock.patch("einsteinpy.plotting.geodesic.static.plt.show")
def test_static_geod_plotter_draws_plot2D(mock_show, dummy_geod):
    geod = dummy_geod

    sgpl = StaticGeodesicPlotter()
    sgpl.plot2D(geod)
    fig = sgpl.show()

    mock_show.assert_called_with()


def test_static_geod_plotter_plot2D_raises_error(dummy_geod):
    geod = dummy_geod

    sgpl = StaticGeodesicPlotter()

    try:
        sgpl.plot2D(geod, coordinates=(-1, 2))

        assert False

    except IndexError:
        assert True


@mock.patch("einsteinpy.plotting.geodesic.static.plt.show")
def test_static_geod_plotter_draws_parametric_plot(mock_show, dummy_geod):
    geod = dummy_geod

    sgpl = StaticGeodesicPlotter()
    sgpl.parametric_plot(geod)
    fig = sgpl.show()

    mock_show.assert_called_with()


def test_static_geod_plotter_parameters(dummy_geod):
    geod = dummy_geod

    sgpl = StaticGeodesicPlotter(bh_colors=("#0F0", "#FAF"), draw_ergosphere=False)
    assert sgpl.bh_colors == ("#0F0", "#FAF")
    assert sgpl.draw_ergosphere is False


def test_static_geod_plotter_ax_warning(dummy_geod):
    geod = dummy_geod

    fig, ax = plt.subplots()

    assert isinstance(ax, Axes)

    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")

        sgpl = StaticGeodesicPlotter(ax=ax)

        assert len(w) == 1  # 1 warning to be shown
        assert issubclass(w[-1].category, PendingDeprecationWarning)


@mock.patch("einsteinpy.plotting.geodesic.static.plt.show")
def test_plot_calls_plt_geod_show(mock_show, dummy_geod):
    geod = dummy_geod

    sgpl = StaticGeodesicPlotter()
    sgpl.plot(geod)
    sgpl.show()

    mock_show.assert_called_with()


@mock.patch("einsteinpy.plotting.geodesic.static.plt.show")
def test_plot_calls_plt_geod_show_view_init(mock_show, dummy_geod):
    geod = dummy_geod

    sgpl = StaticGeodesicPlotter()
    sgpl.plot(geod)
    sgpl.show(azim=-60, elev=30)

    mock_show.assert_called_with()


@mock.patch("einsteinpy.plotting.geodesic.static.plt.savefig")
def test_plot_geod_save_saves_plot(mock_save, dummy_geod):
    geod = dummy_geod

    sgpl = StaticGeodesicPlotter()
    sgpl.plot(geod)
    name = "test_plot.png"
    sgpl.save(name)

    mock_save.assert_called_with(name)
    assert sgpl.ax.name == "3d"


@mock.patch("einsteinpy.plotting.geodesic.static.plt.savefig")
def test_plot_geod_save_saves_parametric_plot(mock_save, dummy_geod):
    geod = dummy_geod

    sgpl = StaticGeodesicPlotter()
    sgpl.parametric_plot(geod)
    name = "Parametric.png"
    sgpl.save()

    mock_save.assert_called_with(name)
    assert sgpl.ax.name != "3d"


@mock.patch("einsteinpy.plotting.geodesic.static.plt.show")
def test_static_geod_plotter_show_clear(mock_show, dummy_geod):
    geod = dummy_geod

    sgpl = StaticGeodesicPlotter(draw_ergosphere=False)
    sgpl.plot(geod)
    sgpl.show()

    mock_show.assert_called_with()
    sgpl.clear()
    assert sgpl.fig.get_axes() == []


def test_animate_creates_ani(dummy_geod):
    geod = dummy_geod

    sgpl = StaticGeodesicPlotter()
    sgpl.animate(geod, interval=10)

    assert sgpl.ani
