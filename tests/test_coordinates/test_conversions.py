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
        t=0.,
        x=10 / np.sqrt(2),
        y=10 / np.sqrt(2),
        z=0.,
        v_x=-190 / np.sqrt(2),
        v_y=210 / np.sqrt(2),
        v_z=200.
    )


@pytest.fixture()
def spherical_coords():
    return SphericalConversion(
        t=0.,
        r=10.,
        theta=np.pi / 2,
        phi=np.pi / 4,
        v_r=10.,
        v_th=-20.,
        v_p=20.
    )


@pytest.fixture()
def bl_coords():
    return BoyerLindquistConversion(
        t=0.,
        r=10.,
        theta=np.pi / 2,
        phi=np.pi / 4,
        v_r=10.,
        v_th=-20.,
        v_p=20.
    )


def strip_velocities(coords):
    vals = coords.values()
    return type(coords)(*vals[:4])


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
    M = 1e24
    a = 0.75
    to_bl_coords = cartesian_coords.convert_bl(M=M, a=a)
    assert_allclose(to_bl_coords, bl_coords.values(), rtol=0.0, atol=1e-6)
    cartesian_coords, bl_coords = (
        strip_velocities(cartesian_coords),
        strip_velocities(bl_coords),
    )
    to_bl_coords = cartesian_coords.convert_bl(M=M, a=a)
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
    M = 1e24
    a = 0.75
    to_bl_coords = spherical_coords.convert_bl(M=M, a=a)
    assert_allclose(to_bl_coords, bl_coords.values(), rtol=0.0, atol=1e-6)
    spherical_coords, bl_coords = (
        strip_velocities(spherical_coords),
        strip_velocities(bl_coords),
    )
    to_bl_coords = spherical_coords.convert_bl(M=M, a=a)
    assert_allclose(to_bl_coords, bl_coords.values(), rtol=0.0, atol=1e-6)


def test_BoyerLindquistToCartesianDifferential(bl_coords, cartesian_coords):
    M = 1e24
    a = 0.75
    to_cartesian_coords = bl_coords.convert_cartesian(M=M, a=a)
    assert_allclose(to_cartesian_coords, cartesian_coords.values(), rtol=0.0, atol=1e-6)
    bl_coords, cartesian_coords = (
        strip_velocities(bl_coords),
        strip_velocities(cartesian_coords),
    )
    to_cartesian_coords = bl_coords.convert_cartesian(M=M, a=a)
    assert_allclose(to_cartesian_coords, cartesian_coords.values(), rtol=0.0, atol=1e-6)


def test_BoyerLindquistToSphericalDifferential(bl_coords, spherical_coords):
    M = 1e24
    a = 0.75
    to_spherical_coords = bl_coords.convert_spherical(M=M, a=a)
    assert_allclose(to_spherical_coords, spherical_coords.values(), rtol=0.0, atol=1e-6)
    bl_coords, spherical_coords = (
        strip_velocities(bl_coords),
        strip_velocities(spherical_coords),
    )
    to_spherical_coords = bl_coords.convert_spherical(M=M, a=a)
    assert_allclose(to_spherical_coords, spherical_coords.values(), rtol=0.0, atol=1e-6)


def test_cycle_BLSphericalDifferential(bl_coords):
    M = 1e24
    a = 0.75
    bl_diff = bl_coords
    sph_diff = bl_diff.convert_spherical(M=M, a=a)
    bl_diff2 = SphericalConversion(*sph_diff).convert_bl(M=M, a=a)
    assert_allclose(bl_diff2, bl_diff.values(), rtol=0.0, atol=1e-6)


def test_cycle_BLCartesianDifferential(bl_coords):
    M = 1e24
    a = 0.75
    bl_diff = bl_coords
    cart_diff = bl_diff.convert_cartesian(M=M, a=a)
    bl_diff2 = CartesianConversion(*cart_diff).convert_bl(M=M, a=a)
    assert_allclose(bl_diff2, bl_diff.values(), rtol=0.0, atol=1e-6)
