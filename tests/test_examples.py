from einsteinpy.bodies import Body
from einsteinpy.coordinates import SphericalDifferential
from einsteinpy.examples import perihelion
from einsteinpy.geodesic import Geodesic
from einsteinpy.metric import Schwarzschild


def test_perihelion_works():
    p = perihelion()
    assert isinstance(p.attractor, Body)
    assert isinstance(p, Geodesic)
    assert p.attractor.name == "BH"


def test_perihelion_works_2():
    p = perihelion()
    assert isinstance(p.metric, Schwarzschild)
    assert isinstance(p.metric.input_coords, SphericalDifferential)
