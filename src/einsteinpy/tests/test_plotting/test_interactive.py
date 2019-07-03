import os
from unittest import mock

import astropy.units as u
import numpy as np
import pytest
from plotly.graph_objs import FigureWidget

from einsteinpy.bodies import Body
from einsteinpy.coordinates import SphericalDifferential
from einsteinpy.geodesic import Geodesic
from einsteinpy.plotting import InteractiveGeodesicPlotter


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


def test_plotlyplotter_has_figurewidget(dummy_data):
    geodesic = dummy_data
    cl = InteractiveGeodesicPlotter()
    assert isinstance(cl.fig, FigureWidget)
    assert cl.attractor_present is False


def test_plot_draws_fig(dummy_data):
    geodesic = dummy_data
    cl = InteractiveGeodesicPlotter()
    cl.plot(geodesic)
    fig = cl.show()
    assert cl.attractor_present
    assert fig


@mock.patch("einsteinpy.plotting.geodesics.interactive.saveplot")
def test_save_saves_plot(mock_save, dummy_data):
    geodesic = dummy_data
    cl = InteractiveGeodesicPlotter()
    cl.plot(geodesic)
    name = "test_plot.png"
    cl.save(name)
    basename, ext = os.path.splitext(name)
    mock_save.assert_called_with(
        cl.fig, image=ext[1:], image_filename=basename, show_link=False
    )
