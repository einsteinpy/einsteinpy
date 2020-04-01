import sys
from io import StringIO

import astropy.units as u
import numpy as np
import pytest
from numpy.testing import assert_allclose

from einsteinpy.coordinates import (
    BoyerLindquistConversion,
    CartesianConversion,
    SphericalConversion,
)


@pytest.fixture()
def cartesian_coords():
    return CartesianConversion(
        10 / np.sqrt(2), 10 / np.sqrt(2), 0, -190 / np.sqrt(2), 210 / np.sqrt(2), 200.0
    )


@pytest.fixture()
def spherical_coords():
    return SphericalConversion(
        10.0, 1.5707963267948966, 0.7853981633974483, 10.0, -20.0, 20.0
    )


@pytest.fixture()
def bl_coords():
    return BoyerLindquistConversion(
        10.0, 1.5707963267948966, 0.7853981633974483, 10.0, -20.0, 20.0, 0.0
    )


def strip_velocities(coords):
    if isinstance(coords, BoyerLindquistConversion):
        vals = coords.values()
        return BoyerLindquistConversion(*vals[:3], a=vals[-1])
    elif isinstance(coords, (CartesianConversion, SphericalConversion)):
        vals = coords.values()
        return type(coords)(*vals[:3])
    return None


def test_CartesianToSpherical(cartesian_coords, spherical_coords):
    to_spherical_coords = cartesian_coords.convert_spherical()
    assert_allclose(to_spherical_coords, spherical_coords.values(), rtol=0.0, atol=1e-6)
    cartesian_coords, spherical_coords = (
        strip_velocities(cartesian_coords),
        strip_velocities(spherical_coords),
    )
    to_spherical_coords = cartesian_coords.convert_spherical()
    assert_allclose(to_spherical_coords, spherical_coords.values(), rtol=0.0, atol=1e-6)


def test_CartesianToBoyerLindquistDifferential(cartesian_coords, bl_coords):
    to_bl_coords = cartesian_coords.convert_bl(bl_coords.a_si)
    assert_allclose(to_bl_coords, bl_coords.values(), rtol=0.0, atol=1e-6)
    cartesian_coords, bl_coords = (
        strip_velocities(cartesian_coords),
        strip_velocities(bl_coords),
    )
    to_bl_coords = cartesian_coords.convert_bl(bl_coords.a_si)
    assert_allclose(to_bl_coords, bl_coords.values(), rtol=0.0, atol=1e-6)


def test_SphericalToCartesianDifferential(spherical_coords, cartesian_coords):
    to_cartesian_coords = spherical_coords.convert_cartesian()
    assert_allclose(to_cartesian_coords, cartesian_coords.values(), rtol=0.0, atol=1e-6)
    cartesian_coords, spherical_coords = (
        strip_velocities(cartesian_coords),
        strip_velocities(spherical_coords),
    )
    to_cartesian_coords = spherical_coords.convert_cartesian()
    assert_allclose(to_cartesian_coords, cartesian_coords.values(), rtol=0.0, atol=1e-6)


def test_SphericalToBoyerLindquistDifferential(spherical_coords, bl_coords):
    to_bl_coords = spherical_coords.convert_bl(bl_coords.a_si)
    assert_allclose(to_bl_coords, bl_coords.values(), rtol=0.0, atol=1e-6)
    spherical_coords, bl_coords = (
        strip_velocities(spherical_coords),
        strip_velocities(bl_coords),
    )
    to_bl_coords = spherical_coords.convert_bl(bl_coords.a_si)
    assert_allclose(to_bl_coords, bl_coords.values(), rtol=0.0, atol=1e-6)


def test_BoyerLindquistToCartesianDifferential(bl_coords, cartesian_coords):
    to_cartesian_coords = bl_coords.convert_cartesian()
    assert_allclose(to_cartesian_coords, cartesian_coords.values(), rtol=0.0, atol=1e-6)
    bl_coords, cartesian_coords = (
        strip_velocities(bl_coords),
        strip_velocities(cartesian_coords),
    )
    to_cartesian_coords = bl_coords.convert_cartesian()
    assert_allclose(to_cartesian_coords, cartesian_coords.values(), rtol=0.0, atol=1e-6)


def test_BoyerLindquistToSphericalDifferential(bl_coords, spherical_coords):
    to_spherical_coords = bl_coords.convert_spherical()
    assert_allclose(to_spherical_coords, spherical_coords.values(), rtol=0.0, atol=1e-6)
    bl_coords, spherical_coords = (
        strip_velocities(bl_coords),
        strip_velocities(spherical_coords),
    )
    to_spherical_coords = bl_coords.convert_spherical()
    assert_allclose(to_spherical_coords, spherical_coords.values(), rtol=0.0, atol=1e-6)


def test_cycle_BLSphericalDifferential(bl_coords):
    bl_diff = bl_coords
    sph_diff = bl_diff.convert_spherical()
    bl_diff2 = SphericalConversion(*sph_diff).convert_bl(bl_diff.a_si)
    assert_allclose(bl_diff2, bl_diff.values(), rtol=0.0, atol=1e-6)


def test_cycle_BLCartesianDifferential(bl_coords):
    bl_diff = bl_coords
    cart_diff = bl_diff.convert_cartesian()
    bl_diff2 = CartesianConversion(*cart_diff).convert_bl(bl_diff.a_si)
    assert_allclose(bl_diff2, bl_diff.values(), rtol=0.0, atol=1e-6)
