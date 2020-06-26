import astropy.units as u
import pytest

from einsteinpy.bodies import (
    Body,
    Earth,
    Jupiter,
    Mars,
    Mercury,
    Moon,
    Neptune,
    Pluto,
    Saturn,
    Sun,
    Uranus,
    Venus,
)
from einsteinpy.coordinates import CartesianDifferential, SphericalDifferential


@pytest.mark.parametrize(
    "obj, parent, R, mass",
    [
        (Mercury, Sun, 2439.7 * u.km, 3.285e23 * u.kg),
        (Venus, Sun, 6051.8 * u.km, 4.867e24 * u.kg),
        (Sun, None, 695510 * u.km, 1.9891e30 * u.kg),
        (Earth, Sun, 6731 * u.km, 5.97219e24 * u.kg),
        (Moon, Earth, 1737.5 * u.km, 7.34767309e22 * u.kg),
        (Mars, Sun, 3389.5 * u.km, 6.39e23 * u.kg),
        (Jupiter, Sun, 69911 * u.km, 1.89813e27 * u.kg),
        (Saturn, Sun, 58232 * u.km, 5.683e26 * u.kg),
        (Uranus, Sun, 25362 * u.km, 8.681e25 * u.kg),
        (Neptune, Sun, 24622 * u.km, 1.024e26 * u.kg),
        (Pluto, Sun, 1183.3 * u.km, 1.309e22 * u.kg),
    ],
)
def test_predefined_bodies_base_properties(obj, parent, R, mass):
    assert obj.parent == parent
    assert obj.R == R
    assert obj.mass == mass


@pytest.mark.parametrize(
    "obj, name, q",
    [
        (Mercury, "Mercury", 0 * u.C),
        (Venus, "Venus", 0 * u.C),
        (Sun, "Sun", 0 * u.C),
        (Earth, "Earth", 0 * u.C),
        (Moon, "Moon", 0 * u.C),
        (Mars, "Mars", 0 * u.C),
        (Jupiter, "Jupiter", 0 * u.C),
        (Saturn, "Saturn", 0 * u.C),
        (Uranus, "Uranus", 0 * u.C),
        (Neptune, "Neptune", 0 * u.C),
        (Pluto, "Pluto", 0 * u.C),
    ],
)
def test_predefined_bodies_extra(obj, name, q):
    assert obj.name == name
    assert obj.q == q


def test_state_vecs():
    parent = Sun
    name = "Earth"
    R = 6731 * u.km
    mass = 5.97219e24 * u.kg
    state_vec1 = CartesianDifferential(
        0 * u.km, 0 * u.km, 0 * u.km, 0 * u.km / u.s, 0 * u.km / u.s, 0 * u.km / u.s
    )
    state_vec2 = SphericalDifferential(
        0 * u.km, 0 * u.rad, 0 * u.rad, 0 * u.km / u.s, 0 * u.rad / u.s, 0 * u.rad / u.s
    )
    a = Body(name=name, mass=mass, R=R, state_vec=state_vec1, parent=parent)
    b = Body(name=name, mass=mass, R=R, state_vec=state_vec2, parent=parent)

    assert isinstance(a.pos_vec, list)
    assert isinstance(a.vel_vec, list)
    assert isinstance(b.pos_vec, list)
    assert isinstance(b.vel_vec, list)


def test_body_str_return():
    body = Body(name="BodyTest", mass=1.989e30 * u.kg, R=30 * u.km)

    assert (
        body.__repr__() == "Body: ( Name: (BodyTest), Mass: (1.989e+30 kg), Charge: (0.0 C)', Radius: (30.0 km), Initial State: (None), Parent Body: (None) )"
    )


def test_body_repr_return():
    body = Body(name="BodyTest", mass=1.989e30 * u.kg, R=30 * u.km)

    assert (
        body.__repr__() == "Body: ( Name: (BodyTest), Mass: (1.989e+30 kg), Charge: (0.0 C)', Radius: (30.0 km), Initial State: (None), Parent Body: (None) )"
    )
