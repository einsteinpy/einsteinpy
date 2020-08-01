import os
from unittest import mock

import numpy as np
import pytest
from plotly.graph_objects import Figure

from einsteinpy.geodesic import Nulllike
from einsteinpy.plotting import InteractiveNullGeodesicPlotter


@pytest.fixture()
def dummy_data():
    pos = np.array([0, 30., np.pi / 2, np.pi / 2])
    vel = np.array([-0.2, 0., 0.002])
    a = 0.9
    end_lambda = 200
    max_steps = 200
    return_cartesian = True

    geod = Nulllike(
        pos,
        vel,
        a,
        end_lambda,
        max_steps,
        return_cartesian
    )

    return geod


def test_interactive_nullgeod_plotter_has_figure(dummy_data):
    geodesic = dummy_data
    cl = InteractiveNullGeodesicPlotter()

    assert isinstance(cl.fig, Figure)


def test_plot_nulgeodl_draws_fig(dummy_data):
    geodesic = dummy_data
    cl = InteractiveNullGeodesicPlotter()
    cl.plot(geodesic)
    fig = cl.show()

    assert fig


@mock.patch("einsteinpy.plotting.geodesics.null.interactive.saveplot")
def test_save_saves_plot_nullgeod(mock_save, dummy_data):
    geodesic = dummy_data
    cl = InteractiveNullGeodesicPlotter()
    cl.plot(geodesic)
    name = "test_plot.png"
    cl.save(name)
    basename, ext = os.path.splitext(name)

    mock_save.assert_called_with(
        cl.fig, image=ext[1:], image_filename=basename, show_link=False
    )
