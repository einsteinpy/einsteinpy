import pytest

from einsteinpy.symbolic import MetricTensor
from einsteinpy.symbolic.predefined import (
    AntiDeSitter,
    AntiDeSitterStatic,
    CMetric,
    DeSitter,
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
        Minkowski(),
        CMetric(),
    ],
)
def test_all_predefined_metrics(metric_instance):
    assert isinstance(metric_instance, MetricTensor)
