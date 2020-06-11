from unittest import mock

import pytest
from astropy import units as u
from matplotlib import pyplot as plt

from einsteinpy.plotting import ShadowPlotter
from einsteinpy.rays import Shadow


@pytest.fixture()
def dummy_data():
    mass = 1 * u.kg
    fov = 30 * u.km
    shadow = Shadow(mass=mass, fov=fov, n_rays=1000)
    return shadow


def test_plotter_has_correct_attributes(dummy_data):
    shadow = dummy_data
    cl = ShadowPlotter(shadow=shadow)
    assert isinstance(cl.shadow, Shadow)
    assert cl.is_intensity_plot


@mock.patch("einsteinpy.plotting.rays.shadow.plt.plot")
def test_plot_calls_plt_plot(mock_patch, dummy_data):
    shadow = dummy_data
    cl = ShadowPlotter(shadow=shadow)
    cl.plot()
    expected = [
        mock.call(cl.shadow.fb1, cl.shadow.intensity, "r"),
        mock.call(cl.shadow.fb2, cl.shadow.intensity, "r"),
    ]
    assert mock_patch.call_args_list == expected
