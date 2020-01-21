import astropy.units as u

from einsteinpy.bodies import Body, Earth, Sun
from einsteinpy.coordinates import CartesianDifferential, SphericalDifferential


def test_sun_base_properties():
    assert Sun.parent is None
    assert Sun.R == 695510 * u.km
    assert Sun.mass == 1.9891e30 * u.kg


def test_sun_extras():
    assert Sun.a == 0
    assert Sun.name == "Sun"
    assert Sun.q == 0 * u.C


def test_earth_base_properties():
    assert Earth.parent == Sun
    assert Earth.R == 6731 * u.km
    assert Earth.mass == 5.97219e24 * u.kg


def test_earth_extras():
    assert Earth.a == 0
    assert Earth.name == "Earth"
    assert Earth.q == 0 * u.C


def test_differentials():
    parent = Sun
    name = "Earth"
    R = 6731 * u.km
    mass = 5.97219e24 * u.kg
    differential1 = CartesianDifferential(
        0 * u.km, 0 * u.km, 0 * u.km, 0 * u.km / u.s, 0 * u.km / u.s, 0 * u.km / u.s
    )
    differential2 = SphericalDifferential(
        0 * u.km, 0 * u.rad, 0 * u.rad, 0 * u.km / u.s, 0 * u.rad / u.s, 0 * u.rad / u.s
    )
    a = Body(name=name, mass=mass, R=R, differential=differential1, parent=parent)
    b = Body(name=name, mass=mass, R=R, differential=differential2, parent=parent)
    assert isinstance(a.pos_vec, list)
    assert isinstance(a.vel_vec, list)
    assert isinstance(b.pos_vec, list)
    assert isinstance(b.vel_vec, list)


def test_body_str_return():
    body = Body(name="BodyTest", mass=1.989e30 * u.kg, a=0.3 * u.m, R=30 * u.km)
    assert (
        body.__str__() == "Body ( name: (BodyTest), "
        "mass: (1.989e+30 kg), radius: (30.0 km), "
        "coordinates: (None), spin factor: (0.3 m), charge: (0.0 C) )"
    )


def test_body_repr_return():
    body = Body(name="BodyTest", mass=1.989e30 * u.kg, a=0.3 * u.m, R=30 * u.km)
    assert (
        body.__repr__() == "'Body ( name: (BodyTest), "
        "mass: (1.989e+30 kg), "
        "radius: (30.0 km), "
        "coordinates: (None), spin factor: (0.3 m), charge: (0.0 C) )'"
    )
