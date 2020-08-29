import os
from unittest import mock

import numpy as np
import pytest
from plotly.graph_objects import Figure

from einsteinpy.geodesic import Geodesic
from einsteinpy.plotting import InteractiveGeodesicPlotter


@pytest.fixture()
def dummy_geod():
    q0 = [2.5, np.pi / 2, 0.]
    p0 = [0., -0.2, -2.]
    a = 0.9
    end_lambda = 1.
    step_size = 0.05
    time_like = False
    return_cartesian = True
    julia = False

    geod = Geodesic(
        position=q0,
        momentum=p0,
        a=a,
        end_lambda=end_lambda,
        step_size=step_size,
        time_like=time_like,
        return_cartesian=return_cartesian,
        julia=julia
    )

    return geod


def test_interactive_geod_plotter_has_figure(dummy_geod):
    geod = dummy_geod
    igpl = InteractiveGeodesicPlotter()

    assert isinstance(igpl.fig, Figure)


def test_interactive_geod_plotter_draws_plot(dummy_geod):
    geod = dummy_geod

    igpl = InteractiveGeodesicPlotter()
    igpl.plot(geod)
    fig = igpl.show()

    assert fig


def test_interactive_geod_plotter_draws_plot2D(dummy_geod):
    geod = dummy_geod

    igpl = InteractiveGeodesicPlotter()
    igpl.plot2D(geod)
    fig = igpl.show()

    assert fig


def test_interactive_geod_plotter_plot2D_raises_error(dummy_geod):
    geod = dummy_geod

    igpl = InteractiveGeodesicPlotter()

    try:
        igpl.plot2D(geod, coordinates=(0, 2))

        assert False

    except IndexError:
        assert True


def test_interactive_geod_plotter_draws_parametric_plot(dummy_geod):
    geod = dummy_geod

    igpl = InteractiveGeodesicPlotter()
    igpl.parametric_plot(geod)
    fig = igpl.show()

    assert fig


def test_interactive_geod_plotter_parameters(dummy_geod):
    geod = dummy_geod

    igpl = InteractiveGeodesicPlotter(bh_colors=("#0F0", "#FAF"), draw_ergosphere=False)
    assert igpl.bh_colors == ("#0F0", "#FAF")
    assert igpl.draw_ergosphere is False
    fig = igpl.show()

    assert fig


def test_interactive_geod_plotter_show_clear(dummy_geod):
    geod = dummy_geod

    igpl = InteractiveGeodesicPlotter(draw_ergosphere=False)
    igpl.plot(geod)

    assert isinstance(igpl.show(), Figure)
    igpl.clear()
    assert igpl.fig.data == ()


@mock.patch("einsteinpy.plotting.geodesic.interactive.saveplot")
def test_interactive_geod_plotter_saves_plot(mock_save, dummy_geod):
    geod = dummy_geod
    igpl = InteractiveGeodesicPlotter()
    igpl.plot(geod)
    name = "test_plot.png"
    igpl.save(name)
    basename, ext = os.path.splitext(name)

    mock_save.assert_called_with(
        igpl.fig, image=ext[1:], image_filename=basename, show_link=False
    )
