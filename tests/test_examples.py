from einsteinpy.examples import perihelion
from einsteinpy.geodesic import TimelikeGeodesic
from einsteinpy.metric import Schwarzschild


def test_perihelion_attr():
    """
    Checks various attributes of perihelion()

    """
    p = perihelion()
    assert isinstance(p, TimelikeGeodesic)
    assert p.state.shape[0] == 8
    assert p.trajectory.shape[1] == 8
    assert isinstance(p.metric, Schwarzschild)
