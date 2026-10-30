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


def test_convert_kwargs_raises_KeyError0(
    cartesian_coords, spherical_coords, bl_coords
):
    def c_s(cartesian_coords):
        return cartesian_coords.convert_spherical()

    def c_b(cartesian_coords):
        return cartesian_coords.convert_bl()

    def s_b(spherical_coords):
        return spherical_coords.convert_bl()

    # This should not raise KeyError
    try:
        c_s(cartesian_coords)
    except KeyError:
        pytest.fail("Unexpected KeyError!")

    # These 2 should raise KeyError
    with pytest.raises(KeyError):
        c_b(cartesian_coords)

    with pytest.raises(KeyError):
        s_b(spherical_coords)


def test_cartesian_conversion_values_without_velocities():
    """CartesianConversion.values() returns 4-tuple when velocities not provided."""
    coords = CartesianConversion(t=0.0, x=1.0, y=2.0, z=3.0)
    vals = coords.values()
    assert len(vals) == 4
    assert vals == (0.0, 1.0, 2.0, 3.0)


def test_base_curvilinear_values_without_velocities():
    """BaseCurvilinearConversion.values() returns 4-tuple when velocities not provided."""
    coords = SphericalConversion(t=0.0, r=10.0, theta=np.pi / 2, phi=np.pi / 4)
    vals = coords.values()
    assert len(vals) == 4
    assert vals == (0.0, 10.0, np.pi / 2, np.pi / 4)


def test_base_curvilinear_values_with_velocities():
    """BaseCurvilinearConversion.values() returns 7-tuple when velocities provided."""
    coords = SphericalConversion(
        t=0.0, r=10.0, theta=np.pi / 2, phi=np.pi / 4,
        v_r=10.0, v_th=-20.0, v_p=20.0
    )
    vals = coords.values()
    assert len(vals) == 7
    assert vals == (0.0, 10.0, np.pi / 2, np.pi / 4, 10.0, -20.0, 20.0)


def test_convert_kwargs_raises_KeyError1(
    cartesian_coords, spherical_coords, bl_coords
):
    def s_c(spherical_coords):
        return spherical_coords.convert_cartesian()

    def b_c(bl_coords):
        return bl_coords.convert_cartesian()

    def b_s(bl_coords):
        return bl_coords.convert_spherical()

    # This should not raise KeyError
    try:
        s_c(spherical_coords)
    except KeyError:
        pytest.fail("Unexpected KeyError!")

    # These 2 should raise KeyError
    with pytest.raises(KeyError):
        b_c(bl_coords)

    with pytest.raises(KeyError):
        b_s(bl_coords)
