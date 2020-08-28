from numpy.testing import assert_allclose

from einsteinpy.examples import perihelion
from einsteinpy.geodesic import Timelike


def test_perihelion_attr():
    """
    Checks various attributes of perihelion()

    """
    p = perihelion()
    L = p.trajectory[1][0, 5]

    assert isinstance(p, Timelike)
    assert p.trajectory[1].shape[1] == 6
    assert_allclose(p.trajectory[1][:, 4], 0, atol=1e-12, rtol=1e-12)
    assert_allclose(p.trajectory[1][:, 5], L, atol=1e-12, rtol=1e-12)
