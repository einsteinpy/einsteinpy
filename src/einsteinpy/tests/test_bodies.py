import astropy.units as u

from einsteinpy.bodies import Earth, Sun


def test_sun_base_properties():
    assert Sun.parent == None
    assert Sun.R == 695510 * u.km
    assert Sun.mass == 1.9891e30 * u.kg


def test_sun_extras():
    assert Sun.a == 0
    assert Sun.id == "Sun"
    assert Sun.q == 0 * u.C


def test_earth_base_properties():
    assert Earth.parent == Sun
    assert Earth.R == 6731 * u.km
    assert Earth.mass == 5.97219e24 * u.kg


def test_earth_extras():
    assert Earth.a == 0
    assert Earth.id == "Earth"
    assert Earth.q == 0 * u.C
