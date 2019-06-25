import sys
from io import StringIO

import astropy.units as u
import numpy as np
import pytest
from numpy.testing import assert_allclose

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
    return coordinates.BoyerLindquistDifferential(
        10.0 * u.km,
        1.5707963267948966 * u.rad,
        0.7853981633974483 * u.rad,
        10.0 * u.km / u.s,
        -20.0 * u.rad / u.s,
        20.0 * u.rad / u.s,
        0.0 * u.km,
    )


@pytest.fixture()
def bl_differential2():
    return coordinates.BoyerLindquistDifferential(
        4 * u.km,
        2 * u.rad,
        1 * u.rad,
        -100 * u.m / u.s,
        -1 * u.rad / u.s,
        10 * u.deg / u.s,
        100 * u.m,
    )


def test_CartesianToSphericalDifferential(
    cartesian_differential, spherical_differential
):
    to_spherical_differential = cartesian_differential.spherical_differential()
    assert_allclose(
        to_spherical_differential.si_values(),
        spherical_differential.si_values(),
        rtol=0.0,
        atol=1e-6,
    )


def test_CartesianToBoyerLindquistDifferential(cartesian_differential, bl_differential):
    to_bl_differential = cartesian_differential.bl_differential(bl_differential.a)
    assert_allclose(
        to_bl_differential.si_values(), bl_differential.si_values(), rtol=0.0, atol=1e-6
    )


def test_SphericalToCartesianDifferential(
    spherical_differential, cartesian_differential
):
    to_cartesian_differential = spherical_differential.cartesian_differential()
    assert_allclose(
        to_cartesian_differential.si_values(),
        cartesian_differential.si_values(),
        rtol=0.0,
        atol=1e-6,
    )


def test_SphericalToBoyerLindquistDifferential(spherical_differential, bl_differential):
    to_bl_differential = spherical_differential.bl_differential(bl_differential.a)
    assert_allclose(
        to_bl_differential.si_values(), bl_differential.si_values(), rtol=0.0, atol=1e-6
    )


def test_BoyerLindquistToCartesianDifferential(bl_differential, cartesian_differential):
    to_cartesian_differential = bl_differential.cartesian_differential()
    assert_allclose(
        to_cartesian_differential.si_values(),
        cartesian_differential.si_values(),
        rtol=0.0,
        atol=1e-6,
    )


def test_BoyerLindquistToSphericalDifferential(bl_differential, spherical_differential):
    to_spherical_differential = bl_differential.spherical_differential()
    assert_allclose(
        to_spherical_differential.si_values(),
        spherical_differential.si_values(),
        rtol=0.0,
        atol=1e-6,
    )


def test_cycle_BLSphericalDifferential(bl_differential2):
    bl_diff = bl_differential2
    sph_diff = bl_diff.spherical_differential()
    bl_diff2 = sph_diff.bl_differential(bl_diff.a)
    assert_allclose(bl_diff2.si_values(), bl_diff.si_values(), rtol=0.0, atol=1e-6)


def test_cycle_BLCartesianDifferential(bl_differential2):
    bl_diff = bl_differential2
    cart_diff = bl_diff.cartesian_differential()
    bl_diff2 = cart_diff.bl_differential(bl_diff.a)
    assert_allclose(bl_diff2.si_values(), bl_diff.si_values(), rtol=0.0, atol=1e-6)


# Tests for object.__repr__ and object.__str__


def test_print_core_objects(
    cartesian_differential, spherical_differential, bl_differential
):
    old_stdout = sys.stdout
    # change stdout
    result = StringIO()
    sys.stdout = result

    print(str(cartesian_differential))
    assert "object at 0x" not in result.getvalue()

    print(str(spherical_differential))
    assert "object at 0x" not in result.getvalue()

    print(str(bl_differential))
    assert "object at 0x" not in result.getvalue()

    # again switch to old stdout
    sys.stdout = old_stdout


# Tests for object.velocities()


def test_velocities(cartesian_differential, spherical_differential, bl_differential):
    def with_numpy_array(differential_obj):
        assert_allclose(
            differential_obj.si_values()[3:],
            differential_obj.velocities(return_np=True),
            1e-10,
            1e-15,
        )

    with_numpy_array(cartesian_differential)
    with_numpy_array(spherical_differential)
    with_numpy_array(bl_differential)


def test_velocities2(cartesian_differential, spherical_differential, bl_differential):
    cd, sd, bd = cartesian_differential, spherical_differential, bl_differential
    assert cd.velocities() == [cd.v_x, cd.v_y, cd.v_z]
    assert sd.velocities() == [sd.v_r, sd.v_t, sd.v_p]
    assert bd.velocities() == [bd.v_r, bd.v_t, bd.v_p]
