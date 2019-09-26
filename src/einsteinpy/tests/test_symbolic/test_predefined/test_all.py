import pytest

from einsteinpy.symbolic import MetricTensor
from einsteinpy.symbolic.predefined import AntiDeSitter, AntiDeSitterStatic, DeSitter


@pytest.mark.parametrize(
    "metric_instance", [AntiDeSitter(), AntiDeSitterStatic(), DeSitter()]
)
def test_all_predefined_metrics(metric_instance):
    assert isinstance(metric_instance, MetricTensor)
