import pytest

from einsteinpy.symbolic import MetricTensor
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
