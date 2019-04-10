import astropy.units as u
import pytest

from einsteinpy import coordinates


@pytest.fixture()
def cartesian():
    return coordinates.Cartesian(20.0 * u.km, 311.0 * u.km, 210.0 * u.km)

@pytest.fixture()
def spherical():
    return coordinates.Spherical(375.79382645275 * u.km, 0.9778376650369 * u.rad, 1.5065760775947 * u.rad)

@pytest.fixture()
def boyerlindquist():
    return coordinates.BoyerLindquist(375.79382645275 * u.km, 0.9778376650369 * u.rad, 1.5065760775947 * u.rad), 0.0 * u.km

def assert_cartesian(cartesian, to_cartesian):
    assert abs(to_cartesian.x - cartesian.x) <= (1e-5) * u.km
    assert abs(to_cartesian.y - cartesian.y) <= (1e-5) * u.km
    assert abs(to_cartesian.z - cartesian.z) <= (1e-5) * u.km

def assert_spherical(spherical, to_spherical):
    assert abs(to_spherical.r - spherical.r) <= (1e-5) * u.km
    assert abs(to_spherical.theta - spherical.theta) <= (1e-5) * u.rad
    assert abs(to_spherical.phi - spherical.phi) <= (1e-5) * u.rad

def assert_boyerlindquist(bl, to_bl):
    assert abs(to_bl.r - bl.r) <= (1e-5) * u.km
    assert abs(to_bl.theta - bl.theta) <= (1e-5) * u.rad
    assert abs(to_bl.phi - bl.phi) <= (1e-5) * u.rad

def test_CartesianToSpherical(cartesian, spherical):
    to_spherical = cartesian.to_spherical()
    assert_spherical(spherical, to_spherical)

def test_CartesianToBoyerLindquist(cartesian, boyerlindquist):
    bl, a = boyerlindquist
    to_bl = cartesian.to_bl(a)
    assert_boyerlindquist(bl, to_bl)

def test_SphericalToCartesian(spherical, cartesian):
    to_cartesian = spherical.to_cartesian()
    assert_cartesian(cartesian, to_cartesian)

def test_SphericalToBoyerLindquist(spherical, boyerlindquist):
    bl, a = boyerlindquist
    to_bl = spherical.to_bl(a)
    assert_boyerlindquist(bl, to_bl)

def test_BoyerLindquistToCartesian(boyerlindquist, cartesian):
    bl, a = boyerlindquist
    to_cartesian = bl.to_cartesian(a)
    assert_cartesian(cartesian, to_cartesian)

def test_BoyerLindquistToSpherical(boyerlindquist, spherical):
    bl, a = boyerlindquist
    to_spherical = bl.to_spherical(a)
    assert_spherical(spherical, to_spherical)
