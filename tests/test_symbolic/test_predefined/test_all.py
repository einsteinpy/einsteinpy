import numpy as np
import pytest
from sympy import Array, symbols

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
