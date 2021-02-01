import numpy as np
import warnings
import pytest
from numpy.testing import assert_allclose

from einsteinpy.geodesic import Geodesic
from einsteinpy.geodesic.utils import _kerr
from einsteinpy.integrators import GeodesicIntegrator
from einsteinpy.integrators.utils import _Z


@pytest.mark.parametrize(
    "order, expected",
    [
        (4, (-1.7024143839193153, 1.3512071919596578)),
        (6, (-1.349343516178727, 1.1746717580893635)),
        (8, (-1.2323658786507714, 1.1161829393253857)),
    ]
)
def test__Z(order, expected):
    assert_allclose(_Z(order), expected, atol=1e-8, rtol=1e-8)


def test_str_repr():
    geodint = GeodesicIntegrator(
        metric=_kerr,
        metric_params=(0.9,),
        q0=[2.15, np.pi / 2, 0.],
        p0=[0., 0., 1.5]
    )

    assert str(geodint) == repr(geodint)


def test_runtime_warning1():
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")

        geod = Geodesic(
            metric="Kerr",
            metric_params=(0.9,),
            position=[2.15, np.pi / 2, 0.],
            momentum=[0., 0., 1.5],
            time_like=True,
            steps=4,
            delta=0.5,
            omega=1.  # Unstable integration
        )

        assert len(w) == 2  # 2 warnings to be shown
        assert issubclass(w[-1].category, RuntimeWarning)


def test_runtime_warning2():
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")

        geod = Geodesic(
            metric="Kerr",
            metric_params=(0.9,),
            position=[2.15, np.pi / 2, 0.],
            momentum=[0., 0., 1.5],
            time_like=True,
            steps=4,
            delta=0.5,
            omega=0.01  # Stable integration
        )

        assert len(w) == 0


def test_suppress_runtime_warning():
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")

        geod = Geodesic(
            metric="Kerr",
            metric_params=(0.9,),
            position=[2.15, np.pi / 2, 0.],
            momentum=[0., 0., 1.5],
            time_like=True,
            steps=4,
            delta=0.5,
            omega=1.,  # Unstable integration
            suppress_warnings=True
        )

        assert len(w) == 0


def test_rtol_atol_runtime_warning():
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")

        geod = Geodesic(
            metric="Kerr",
            metric_params=(0.9,),
            position=[2.15, np.pi / 2, 0.],
            momentum=[0., 0., 1.5],
            time_like=True,
            steps=4,
            delta=0.5,
            omega=1.,  # Unstable integration
            rtol=1.,
            atol=1.
        )

        assert len(w) == 0  # 1 warning to be shown


def test_order_NotImplementedError():

    with pytest.raises(NotImplementedError):
        geod = Geodesic(
            metric="Kerr",
            metric_params=(0.9,),
            position=[2.15, np.pi / 2, 0.],
            momentum=[0., 0., 1.5],
            time_like=True,
            steps=4,
            delta=0.5,
            order=5
        )


@pytest.mark.parametrize("order", [4, 6, 8])
def test_higher_order_traits(order):
    geod = Geodesic(
        metric="Kerr",
        metric_params=(0.5,),
        position=[4., np.pi / 2, 0.],
        momentum=[0., 0., 2.],
        time_like=False,
        steps=10,
        delta=0.5,
        order=order,
        return_cartesian=False,
        suppress_warnings=True
    )

    L = geod.momentum[-1]
    theta = geod.position[2]

    assert_allclose(geod.trajectory[1][:, -1], L, atol=1e-4, rtol=1e-4)
    assert_allclose(geod.trajectory[1][:, 2], theta, atol=1e-6, rtol=1e-6)
