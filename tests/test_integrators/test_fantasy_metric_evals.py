import numpy as np
import pytest
from numpy.testing import assert_allclose

from einsteinpy.geodesic import Timelike
from einsteinpy.geodesic.utils import _kerr, _kerrnewman, _sch
from einsteinpy.integrators import GeodesicIntegrator
from einsteinpy.utils.dual import _diff_g, _jacobian_g

# Fixed initial conditions, used for both the reference values and the
# evaluation counts below
_CASES = {
    "Schwarzschild": dict(
        metric="Schwarzschild",
        metric_params=(),
        position=[40.0, np.pi / 2, 0.0],
        momentum=[0.0, 0.0, 3.83405],
    ),
    "Kerr": dict(
        metric="Kerr",
        metric_params=(0.6,),
        position=[10.0, np.pi / 3, 0.0],
        momentum=[0.0, 0.3, 3.2],
    ),
    "KerrNewman": dict(
        metric="KerrNewman",
        metric_params=(0.5, 0.3),
        position=[20.0, np.pi / 3, 0.0],
        momentum=[0.0, 0.5, 4.0],
    ),
}

# State (q1, p1) after 30 steps (delta = 0.5, omega = 1.0), generated with
# einsteinpy at commit 677dfef, before the metric-evaluation changes
_REFERENCE = {
    ("Schwarzschild", 2): [
        15.460524333750278,
        39.95358148436564,
        1.5707963267948966,
        0.03597202434909404,
        -0.9791466489146937,
        -0.006527730797612501,
        8.445085661208105e-18,
        3.8340499987707997,
    ],
    ("Schwarzschild", 4): [
        15.460524623600461,
        39.953574181739214,
        1.5707963267948966,
        0.035972058809511945,
        -0.9791466044305606,
        -0.006516932130697831,
        8.445095047264425e-18,
        3.8340500004448295,
    ],
    ("Kerr", 2): [
        17.868030298989126,
        10.075643979429664,
        1.1717163379329811,
        0.6232961157742088,
        -0.9576672504377013,
        0.012351878449423686,
        1.2980462902844083,
        3.2000039974891057,
    ],
    ("Kerr", 4): [
        17.868024208216674,
        10.07565242109223,
        1.1716992353885831,
        0.6232772298923155,
        -0.9576668211494861,
        0.012330278668315807,
        1.29807521001791,
        3.200000736338456,
    ],
    ("KerrNewman", 2): [
        16.230292581203756,
        19.986182396213522,
        1.074394696925966,
        0.1992451775322245,
        -0.9745563268889069,
        -0.0020459025746354703,
        0.9432438253220594,
        3.999999870436531,
    ],
    ("KerrNewman", 4): [
        16.23029295910086,
        19.986178670666114,
        1.0743969320646998,
        0.19924450049672088,
        -0.9745562940071433,
        -0.0020423261802891566,
        0.9432447510414752,
        3.99999999065312,
    ],
}


@pytest.mark.parametrize("name, order", sorted(_REFERENCE))
def test_trajectory_matches_reference(name, order):
    geod = Timelike(
        return_cartesian=False,
        steps=30,
        delta=0.5,
        order=order,
        omega=1.0,
        suppress_warnings=True,
        **_CASES[name],
    )
    final = geod.trajectory[1][-1]

    assert_allclose(final, _REFERENCE[(name, order)], rtol=1e-12, atol=1e-14)


class _CountingMetric:
    def __init__(self, g):
        self.g = g
        self.calls = 0
        self.__name__ = g.__name__

    def __call__(self, x_vec, *params):
        self.calls += 1
        return self.g(x_vec, *params)


@pytest.mark.parametrize(
    "g, g_prms, q0, p0",
    [
        (_sch, (0.0,), [0.0, 40.0, np.pi / 2, 0.0], [-0.98, 0.0, 0.0, 3.83405]),
        (_kerr, (0.6,), [0.0, 10.0, np.pi / 3, 0.0], [-0.95, 0.0, 0.3, 3.2]),
        (_kerrnewman, (0.5, 0.3), [0.0, 20.0, np.pi / 3, 0.0], [-0.97, 0.0, 0.5, 4.0]),
    ],
)
@pytest.mark.parametrize("order, evals_per_step", [(2, 20), (4, 60)])
def test_metric_evaluations_per_step(g, g_prms, q0, p0, order, evals_per_step):
    # A second-order step is made of four Hamiltonian flows. Each flow needs
    # the metric once on floats and once per coordinate on dual numbers, i.e.
    # 5 evaluations, so 20 per step. Higher orders compose 3 lower-order steps.
    metric = _CountingMetric(g)
    geodint = GeodesicIntegrator(
        metric=metric,
        metric_params=g_prms,
        q0=np.array(q0),
        p0=np.array(p0),
        steps=3,
        order=order,
        suppress_warnings=True,
    )

    for _ in range(3):
        geodint.step()

    assert metric.calls == 3 * evals_per_step


@pytest.mark.parametrize(
    "g, g_prms",
    [(_sch, (0.0,)), (_kerr, (0.9,)), (_kerrnewman, (0.5, 0.3))],
)
@pytest.mark.parametrize("wrt", [0, 1, 2, 3])
def test_jacobian_g_matches_diff_g(g, g_prms, wrt):
    coords = [0.3, 6.5, 1.1, 0.7]
    J = _jacobian_g(g, g_prms, coords, wrt)

    expected = np.array(
        [[_diff_g(g, g_prms, coords, (i, j), wrt) for j in range(4)] for i in range(4)]
    )

    assert_allclose(J, expected, rtol=0, atol=0)


def test_jacobian_g_wrt_out_of_range():
    with pytest.raises(ValueError):
        _jacobian_g(_kerr, (0.9,), [0.0, 6.5, 1.1, 0.7], 4)
