import numpy as np
import pytest
from sympy import Array, Rational, diff, simplify, symbols, tanh

from einsteinpy.symbolic import MetricTensor, RicciTensor, simplify_sympy_array
from einsteinpy.symbolic.predefined import (
    AlcubierreWarp,
    AntiDeSitter,
    AntiDeSitterStatic,
    BarriolaVilekin,
    BertottiKasner,
    BesselGravitationalWave,
    CMetric,
    Davidson,
    DeSitter,
    Ernst,
    FriedmannLemaitreRobertsonWalker,
    Godel,
    JanisNewmanWinicour,
    Kerr,
    KerrNewman,
    KerrNewmanDeSitter,
    Minkowski,
    MinkowskiCartesian,
    MinkowskiPolar,
    ReissnerNordstorm,
    Schwarzschild,
)


@pytest.mark.parametrize(
    "metric_instance",
    [
        AlcubierreWarp(),
        AntiDeSitter(),
        AntiDeSitterStatic(),
        DeSitter(),
        Schwarzschild(),
        Schwarzschild(c=1, sch=2),
        Kerr(),
        KerrNewman(),
        Minkowski(),
        MinkowskiCartesian(),
        MinkowskiPolar(),
        CMetric(),
        Davidson(),
        Godel(),
        BesselGravitationalWave(),
        BarriolaVilekin(),
        BertottiKasner(),
        Ernst(),
        JanisNewmanWinicour(),
        ReissnerNordstorm(),
        FriedmannLemaitreRobertsonWalker(),
        KerrNewmanDeSitter(),
    ],
)
def test_all_predefined_metrics(metric_instance):
    assert isinstance(metric_instance, MetricTensor)


@pytest.mark.parametrize(
    "m1, m2",
    [
        (Schwarzschild(), Kerr(a=0)),  # Schwarzschild is a special case of Kerr
        (Kerr(), KerrNewman(Q=0)),  # Kerr is a special case of Kerr-Newman
        (
            ReissnerNordstorm(),
            KerrNewman(a=0),
        ),  # Reissner-Nordstorm is a special case of Kerr-Newman
        (
            Schwarzschild(),
            ReissnerNordstorm(Q=0),
        ),  # Schwarzschild is a special case of Reissner-Nordstorm
    ],
)
def test_check_two_metrics_are_equal(m1, m2):
    zero_arr = Array(np.zeros(shape=m1.tensor().shape, dtype=int))
    assert simplify_sympy_array(m1.tensor() - m2.tensor()) == zero_arr


def test_Minkowski_equality():
    # Minkowski and MinkowskiCartesian are same
    assert simplify_sympy_array(MinkowskiCartesian().tensor()) == simplify_sympy_array(
        Minkowski().tensor()
    )


def test_BesselGravitationalWave_is_vacuum():
    # Kramer (1999) is a vacuum solution, so the Ricci tensor must vanish
    metric = BesselGravitationalWave()
    zero_arr = Array(np.zeros(shape=metric.tensor().shape, dtype=int))
    assert RicciTensor.from_metric(metric).tensor() == zero_arr
    # The wave depends on time through U = C J0(rho) cos(t)
    assert symbols("t") in metric.tensor()[3, 3].free_symbols


def test_AlcubierreWarp_shape_function():
    # Alcubierre (1994) uses r_s = sqrt((x - x_s)^2 + y^2 + z^2), so the
    # bubble is symmetric under rotations about the x axis
    t, x, y, z = symbols("t x y z")
    sigma, R, v = 2, 1, Rational(1, 2)
    metric = AlcubierreWarp(x_s=0, sigma=sigma, R=R, v=v).tensor()
    for i, j in [(0, 0), (0, 1)]:
        assert simplify(y * diff(metric[i, j], z) - z * diff(metric[i, j], y)) == 0
    # On the bubble wall r_s = R, f = tanh(2 sigma R) / (2 tanh(sigma R))
    f_wall = tanh(2 * sigma * R) / (2 * tanh(sigma * R))
    point = {x: 0, y: Rational(3, 5), z: Rational(4, 5)}
    assert simplify(metric[0, 1].subs(point) + v * f_wall) == 0
