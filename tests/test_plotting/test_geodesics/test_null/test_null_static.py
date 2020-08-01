from unittest import mock

import numpy as np
import pytest
from matplotlib.axes import Axes

from einsteinpy.geodesic import Nulllike
from einsteinpy.plotting import StaticNullGeodesicPlotter


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


def test_static_nullgeod_plotter_has_axes(dummy_data):
    geodesic = dummy_data
    cl = StaticNullGeodesicPlotter()

    assert isinstance(cl.ax, Axes)


@mock.patch("einsteinpy.plotting.geodesics.null.static.plt.show")
def test_plot_calls_plt_nullgeod_show(mock_show, dummy_data):
    geodesic = dummy_data
    cl = StaticNullGeodesicPlotter()
    cl.plot(geodesic)
    cl.show()

    mock_show.assert_called_with()


@mock.patch("einsteinpy.plotting.geodesics.null.static.plt.savefig")
def test_plot_nullgeod_save_saves_plot(mock_save, dummy_data):
    geodesic = dummy_data
    cl = StaticNullGeodesicPlotter()
    cl.plot(geodesic)
    name = "test_plot.png"
    cl.save(name)

    mock_save.assert_called_with(name)
    assert cl.ax.name == "3d"
