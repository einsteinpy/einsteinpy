import numpy as np
import pytest
from sympy import Array

from einsteinpy.symbolic import MetricTensor, simplify_sympy_array
from einsteinpy.symbolic.predefined import (
    AntiDeSitter,
    AntiDeSitterStatic,
    BarriolaVilekin,
    BertottiKasner,
    BesselGravitationalWave,
    CMetric,
    Davidson,
    DeSitter,
    Ernst,
    Godel,
    JanisNewmanWinicour,
    Kerr,
    KerrNewman,
    Minkowski,
    MinkowskiCartesian,
    MinkowskiPolar,
    ReissnerNordstorm,
    Schwarzschild,
)


@pytest.mark.parametrize(
    "metric_instance",
    [
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
