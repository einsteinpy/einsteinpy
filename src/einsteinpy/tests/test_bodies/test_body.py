from unittest import mock

import astropy.units as u
import numpy as np
import pytest

from einsteinpy.bodies import Body
from einsteinpy.plotting import ScatterGeodesicPlotter


@pytest.fixture()
def dummy_data():
    tempbody = Body(
        mass=4e24 * u.kg,
        pos_vec=[306 * u.m, np.pi / 2 * u.rad, np.pi / 2 * u.rad],
        vel_vec=[0 * u.m / u.s, 0 * u.rad / u.s, 951.0 * u.rad / u.s],
        desc="This is a temporary body under testing!",
    )
    t = 0 * u.s
    start_lambda = 0.0
    end_lambda = 0.002
    step_size = 0.5e-6
    return (
        tempbody.pos_vec,
        tempbody.vel_vec,
        t,
        tempbody.mass,
        start_lambda,
        end_lambda,
        step_size,
    )


def test_plot_attractor_is_called_only_once(dummy_data):
    r, v, _, m, _, el, ss = dummy_data
    cl = ScatterGeodesicPlotter(m)
    assert not cl._attractor_present
    cl.plot(r, v, el, ss)
    assert cl._attractor_present


@mock.patch(
    "einsteinpy.plotting.geodesics_static.ScatterGeodesicPlotter._plot_attractor"
)
def test_plot_calls_plot_attractor(mock_plot_attractor):
    tempbody = Body(
        mass=4e24 * u.kg,
        pos_vec=[306 * u.m, np.pi / 2 * u.rad, np.pi / 2 * u.rad],
        vel_vec=[0 * u.m / u.s, 0 * u.rad / u.s, 951.0 * u.rad / u.s],
        desc="This is a temporary body under testing!",
    )
    el = 0.002
    ss = 0.5e-6
    cl = ScatterGeodesicPlotter(tempbody.mass)
    cl.plot(tempbody.pos_vec, tempbody.vel_vec, el, ss)
    mock_plot_attractor.assert_called_with()


@mock.patch("einsteinpy.plotting.geodesics_static.plt.show")
def test_plot_show_shows_plot(mock_show):
    tempbody = Body(
        mass=4e24 * u.kg,
        pos_vec=[306 * u.m, np.pi / 2 * u.rad, np.pi / 2 * u.rad],
        vel_vec=[0 * u.m / u.s, 0 * u.rad / u.s, 951.0 * u.rad / u.s],
        desc="This is a temporary body under testing!",
    )
    el = 0.002
    ss = 0.5e-6
    cl = ScatterGeodesicPlotter(tempbody.mass)
    cl.plot(tempbody.pos_vec, tempbody.vel_vec, el, ss)
    cl.show()
    mock_show.assert_called_with()


@mock.patch("einsteinpy.plotting.geodesics_static.plt.savefig")
def test_plot_save_saves_plot(mock_save):
    tempbody = Body(
        mass=4e24 * u.kg,
        pos_vec=[306 * u.m, np.pi / 2 * u.rad, np.pi / 2 * u.rad],
        vel_vec=[0 * u.m / u.s, 0 * u.rad / u.s, 951.0 * u.rad / u.s],
        desc="This is a temporary body under testing!",
    )
    el = 0.002
    ss = 0.5e-6
    cl = ScatterGeodesicPlotter(tempbody.mass)
    cl.plot(tempbody.pos_vec, tempbody.vel_vec, el, ss)
    name = "test_plot.png"
    cl.save(name)
    mock_save.assert_called_with(name)
