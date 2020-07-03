import os
from unittest import mock

import astropy.units as u
import numpy as np
import pytest
from plotly.graph_objects import Figure

from einsteinpy.coordinates import SphericalDifferential
from einsteinpy.metric import Schwarzschild
from einsteinpy.geodesic import TimelikeGeodesic
from einsteinpy.plotting import InteractiveGeodesicPlotter


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


def test_interactive_plotter_has_figure(dummy_data):
    geodesic = dummy_data
    cl = InteractiveGeodesicPlotter()
    assert isinstance(cl.fig, Figure)
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


def test_plot_calls_draw_attractor_Manualscale(dummy_data):
    geodesic = dummy_data
    cl = InteractiveGeodesicPlotter(attractor_radius_scale=1500)
    cl.plot(geodesic)
    assert cl.attractor_present
    assert cl.attractor_radius_scale == 1500


def test_plot_calls_draw_attractor_AutoScale(dummy_data):
    geodesic = dummy_data
    cl = InteractiveGeodesicPlotter()
    cl.plot(geodesic)
    assert cl.attractor_present
