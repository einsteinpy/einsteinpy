import numpy as np
import pytest
from sympy import Array

from einsteinpy.symbolic import MetricTensor, simplify_sympy_array
from einsteinpy.symbolic.predefined import (
    AntiDeSitter,
    AntiDeSitterStatic,
    CMetric,
    Davidson,
    DeSitter,
    Godel,
    Kerr,
    KerrNewman,
    Minkowski,
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
        CMetric(),
        Davidson(),
        Godel(),
    ],
)
def test_all_predefined_metrics(metric_instance):
    assert isinstance(metric_instance, MetricTensor)


@pytest.mark.parametrize(
    "m1, m2",
    [
        (Schwarzschild(), Kerr(a=0)),  # Schwarzschild is a special case of Kerr
        (Kerr(), KerrNewman(Q=0)),  # Kerr is a special case of Kerr-Newman
    ],
)
def test_check_two_metrics_are_equal(m1, m2):
    zero_arr = Array(np.zeros(shape=m1.tensor().shape, dtype=int))
    assert simplify_sympy_array(m1.tensor() - m2.tensor()) == zero_arr
