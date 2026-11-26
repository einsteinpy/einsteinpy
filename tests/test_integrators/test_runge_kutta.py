import warnings

import numpy as np
import pytest

from einsteinpy import integrators


def test_RK45_RuntimeWarning():
    def fun(t, y):
        return 2 * y

    cl = integrators.RK45(fun, 3.0, np.array([3, 2, 1.0]), 8, 2.0)
    with warnings.catch_warnings(record=True) as w:
        for i in range(10):
            cl.step()
        assert len(w) == 1 and issubclass(w[-1].category, RuntimeWarning)


@pytest.mark.parametrize(
    "t0, y0, t_bound, stepsize, y_end",
    [
        (0.125, 1.0, 1.125, 0.25, 2.0),
        (1.125, 1.0, 0.125, -0.25, 0.0),
    ],
)
def test_RK4naive_stops_at_t_bound(t0, y0, t_bound, stepsize, y_end):
    # y' = 1, so y - t is constant; all times here are exact in binary
    cl = integrators.RK4naive(lambda t, y: 1.0, t0, y0, t_bound, stepsize)
    for _ in range(4):
        cl.step()
    assert (cl.t, cl.y) == (t_bound, y_end)
    with pytest.warns(RuntimeWarning):
        cl.step()
    assert (cl.t, cl.y) == (t_bound, y_end)


def test_RK4naive_does_not_step_when_t0_is_t_bound():
    cl = integrators.RK4naive(lambda t, y: 1.0, 2.0, 3.0, 2.0, 0.5)
    with pytest.warns(RuntimeWarning):
        cl.step()
    assert (cl.t, cl.y) == (2.0, 3.0)


def test_RK4naive_exponential():
    cl = integrators.RK4naive(lambda t, y: y, 0.0, 1.0, 1.0, 0.1)
    for _ in range(10):
        cl.step()
    assert np.isclose(cl.t, 1.0)
    assert abs(cl.y - np.e) < 3e-6


@pytest.mark.parametrize(
    "t0, t_bound, stepsize",
    [
        (0.0, 1.0, 0.375),
        (1.0, 0.0, -0.375),
        (1.0, 0.0, 0.375),
        (0.0, 1.0, -0.375),
    ],
)
def test_RK4naive_last_step_lands_on_t_bound(t0, t_bound, stepsize):
    # y' = 1 is integrated exactly, so y - t stays at y0 - t0
    cl = integrators.RK4naive(lambda t, y: 1.0, t0, 5.0, t_bound, stepsize)
    for _ in range(3):
        cl.step()
    assert cl.t == t_bound
    assert cl.y == 5.0 + (t_bound - t0)
    with pytest.warns(RuntimeWarning):
        cl.step()
    assert cl.t == t_bound


def test_RK4naive_stages_stay_inside_interval():
    # stage times must not go past t_bound
    def fun(t, y):
        assert 0.0 <= t <= 1.0
        return np.cos(t) * y

    cl = integrators.RK4naive(fun, 0.0, 1.0, 1.0, 0.3)
    for _ in range(4):
        cl.step()
    assert cl.t == 1.0
    assert abs(cl.y - np.exp(np.sin(1.0))) < 1e-4
