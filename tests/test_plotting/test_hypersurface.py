from unittest import mock

import pytest
from astropy import units as u

from einsteinpy.hypersurface import SchwarzschildEmbedding
from einsteinpy.plotting import HypersurfacePlotter


@pytest.fixture()
def dummy_data():
    surface_obj = SchwarzschildEmbedding(5.927e23 * u.kg)
    return surface_obj


@mock.patch("einsteinpy.plotting.hypersurface.core.plt.show")
def test_plot_calls_plt_show(mock_show, dummy_data):
    surface = dummy_data
    cl = HypersurfacePlotter(surface)
    cl.plot()
    cl.show()
    mock_show.assert_called_with()
    assert cl.alpha == 100
    assert cl.plot_type == "wireframe"


@mock.patch("einsteinpy.plotting.hypersurface.core.HypersurfacePlotter.show")
@mock.patch("einsteinpy.plotting.hypersurface.core.HypersurfacePlotter.plot")
def test_plot_works_with_different_plot_type(mock_plot, mock_show, dummy_data):
    surface = dummy_data
    cl = HypersurfacePlotter(surface, plot_type="surface")
    cl.plot()
    cl.show()
    assert cl.plot_type == "surface"
