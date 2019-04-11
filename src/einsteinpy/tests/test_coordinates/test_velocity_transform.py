import astropy.units as u
import numpy as np
import pytest

from einsteinpy import coordinates


@pytest.fixture()
def cartesian_differential():
    return coordinates.CartesianDifferential(
        10 / np.sqrt(2) * u.km,
        10 / np.sqrt(2) * u.km,
        0 * u.km,
        -190 / np.sqrt(2) * u.km / u.s,
        210 / np.sqrt(2) * u.km / u.s,
        200.0 * u.km / u.s,
    )


@pytest.fixture()
def spherical_differential():
    return coordinates.SphericalDifferential(
        10.0 * u.km,
        1.5707963267948966 * u.rad,
        0.7853981633974483 * u.rad,
        10.0 * u.km / u.s,
        -20.0 * u.rad / u.s,
        20.0 * u.rad / u.s,
    )


@pytest.fixture()
def bl_differential():
    return (
        coordinates.BoyerLindquistDifferential(
            10.0 * u.km,
            1.5707963267948966 * u.rad,
            0.7853981633974483 * u.rad,
            10.0 * u.km / u.s,
            -20.0 * u.rad / u.s,
            20.0 * u.rad / u.s,
        ),
        0.0 * u.km,
    )


def assert_cartesian_pos(cartesian_differential, to_cartesian_differential):
    assert abs(to_cartesian_differential.y - cartesian_differential.y) <= (1e-5) * u.km
    assert abs(to_cartesian_differential.z - cartesian_differential.z) <= (1e-5) * u.km
    assert (
        abs(to_cartesian_differential.v_x - cartesian_differential.v_x)
        <= (1e-5) * u.km / u.s
    )


def assert_cartesian_differential(cartesian_differential, to_cartesian_differential):
    assert abs(to_cartesian_differential.x - cartesian_differential.x) <= (1e-5) * u.km

    assert (
        abs(to_cartesian_differential.v_y - cartesian_differential.v_y)
        <= (1e-5) * u.km / u.s
    )
    assert (
        abs(to_cartesian_differential.v_z - cartesian_differential.v_z)
        <= (1e-5) * u.km / u.s
    )


def assert_spherical_pos(spherical_differential, to_spherical_differential):
    assert abs(to_spherical_differential.r - spherical_differential.r) <= (1e-5) * u.km
    assert (
        abs(to_spherical_differential.theta - spherical_differential.theta)
        <= (1e-5) * u.rad
    )
    assert (
        abs(to_spherical_differential.phi - spherical_differential.phi)
        <= (1e-5) * u.rad
    )


def assert_spherical_differential(spherical_differential, to_spherical_differential):
    assert (
        abs(to_spherical_differential.v_r - spherical_differential.v_r)
        <= (1e-5) * u.km / u.s
    )
    assert (
        abs(to_spherical_differential.v_t - spherical_differential.v_t)
        <= (1e-5) * u.rad / u.s
    )
    assert (
        abs(to_spherical_differential.v_p - spherical_differential.v_p)
        <= (1e-5) * u.rad / u.s
    )


def assert_boyerlindquist_pos(bl_differential, to_bl_differential):
    assert abs(to_bl_differential.r - bl_differential.r) <= (1e-5) * u.km
    assert abs(to_bl_differential.theta - bl_differential.theta) <= (1e-5) * u.rad
    assert abs(to_bl_differential.phi - bl_differential.phi) <= (1e-5) * u.rad


def assert_bl_differential(bl_differential, to_bl_differential):

    assert abs(to_bl_differential.v_r - bl_differential.v_r) <= (1e-5) * u.km / u.s
    assert abs(to_bl_differential.v_t - bl_differential.v_t) <= (1e-5) * u.rad / u.s
    assert abs(to_bl_differential.v_p - bl_differential.v_p) <= (1e-5) * u.rad / u.s


def test_CartesianToSphericalDifferential(
    cartesian_differential, spherical_differential
):
    to_spherical_differential = cartesian_differential.spherical_differential()
    assert_spherical_pos(spherical_differential, to_spherical_differential)
    assert_spherical_differential(spherical_differential, to_spherical_differential)


def test_CartesianToBoyerLindquistDifferential(cartesian_differential, bl_differential):
    bl_differential, a = bl_differential
    to_bl_differential = cartesian_differential.bl_differential(a)
    assert_boyerlindquist_pos(bl_differential, to_bl_differential)
    assert_bl_differential(bl_differential, to_bl_differential)


def test_SphericalToCartesianDifferential(
    spherical_differential, cartesian_differential
):
    to_cartesian_differential = spherical_differential.cartesian_differential()
    assert_cartesian_pos(cartesian_differential, to_cartesian_differential)
    assert_cartesian_differential(cartesian_differential, to_cartesian_differential)


def test_SphericalToBoyerLindquistDifferential(spherical_differential, bl_differential):
    bl_differential, a = bl_differential
    to_bl_differential = spherical_differential.bl_differential(a)
    assert_boyerlindquist_pos(bl_differential, to_bl_differential)
    assert_bl_differential(bl_differential, to_bl_differential)


def test_BoyerLindquistToCartesianDifferential(bl_differential, cartesian_differential):
    bl_differential, a = bl_differential
    to_cartesian_differential = bl_differential.cartesian_differential(a)
    assert_cartesian_pos(cartesian_differential, to_cartesian_differential)
    assert_cartesian_differential(cartesian_differential, to_cartesian_differential)


def test_BoyerLindquistToSphericalDifferential(bl_differential, spherical_differential):
    bl_differential, a = bl_differential
    to_spherical_differential = bl_differential.spherical_differential(a)
    assert_spherical_pos(spherical_differential, to_spherical_differential)
    assert_spherical_differential(spherical_differential, to_spherical_differential)
