from einsteinpy.examples import perihelion
from einsteinpy.geodesic import Geodesic
from einsteinpy.metric import Schwarzschild


def test_perihelion_attr():
    """
    Checks various attributes of perihelion()

    """
    p = perihelion()
    assert isinstance(p, Geodesic)
    assert p.init_vec.shape[0] == 8
    assert p.trajectory.shape[1] == 8
    assert isinstance(p.metric, Schwarzschild)
