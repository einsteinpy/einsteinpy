from unittest import mock

import numpy as np
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


def test_shadow_smoothen(dummy_data):
    """Shadow.smoothen() interpolates intensity for smoother plots."""
    shadow = dummy_data
    fb1_before = np.asarray(shadow.fb1)

    shadow.smoothen(points=100)

    assert len(shadow.fb1) == 100
    assert len(shadow.intensity) == 100
    assert shadow.fb2.shape == shadow.fb1.shape
    assert shadow.fb1.min() == pytest.approx(float(fb1_before.min()))
    assert shadow.fb1.max() == pytest.approx(float(fb1_before.max()))
