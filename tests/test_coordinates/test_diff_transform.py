import sys
from io import StringIO
from unittest import mock

import astropy.units as u
import numpy as np
import pytest
from numpy.testing import assert_allclose

from einsteinpy import coordinates
from einsteinpy.metric import Schwarzschild, Kerr


@pytest.fixture()
def cartesian_differential():
    return coordinates.CartesianDifferential(
        1e3 * u.s,
        10e3 / np.sqrt(2) * u.m,
        10e3 / np.sqrt(2) * u.m,
        0 * u.m,
        -190e3 / np.sqrt(2) * u.m / u.s,
        210e3 / np.sqrt(2) * u.m / u.s,
        200.0e3 * u.m / u.s
    )


@pytest.fixture()
def spherical_differential():
    return coordinates.SphericalDifferential(
        1e3 * u.s,
        10e3 * u.m,
        1.5707963267948966 * u.rad,
        0.7853981633974483 * u.rad,
        10e3 * u.m / u.s,
        -20.0 * u.rad / u.s,
        20.0 * u.rad / u.s
    )


@pytest.fixture()
def bl_differential():
    return coordinates.BoyerLindquistDifferential(
        1e3 * u.s,
        10e3 * u.m,
        1.5707963267948966 * u.rad,
        0.7853981633974483 * u.rad,
        10e3 * u.m / u.s,
        -20.0 * u.rad / u.s,
        20.0 * u.rad / u.s
    )


@pytest.fixture()
def bl_differential2():
    return coordinates.BoyerLindquistDifferential(
        1e3 * u.s,
        4e3 * u.m,
        2 * u.rad,
        1 * u.rad,
        -100 * u.m / u.s,
        -1 * u.rad / u.s,
        0.17453292519943295 * u.rad / u.s
    )


def test_CartesianToSphericalDifferential(
    cartesian_differential, spherical_differential
):
    to_spherical_differential = cartesian_differential.spherical_differential()
    assert_allclose(
        to_spherical_differential.values(),
        spherical_differential.values(),
        rtol=0.0,
        atol=1e-6,
    )


def test_CartesianToBoyerLindquistDifferential(cartesian_differential, bl_differential):
    M = 1e24
    a = 0.75
    to_bl_differential = cartesian_differential.bl_differential(M=M, a=a)
    assert_allclose(
        to_bl_differential.values(), bl_differential.values(), rtol=0.0, atol=1e-6
    )


def test_SphericalToCartesianDifferential(
    spherical_differential, cartesian_differential
):
    to_cartesian_differential = spherical_differential.cartesian_differential()
    assert_allclose(
        to_cartesian_differential.values(),
        cartesian_differential.values(),
        rtol=0.0,
        atol=1e-6,
    )


def test_SphericalToBoyerLindquistDifferential(spherical_differential, bl_differential):
    M = 1e24
    a = 0.75
    to_bl_differential = spherical_differential.bl_differential(M=M, a=a)
    assert_allclose(
        to_bl_differential.values(), bl_differential.values(), rtol=0.0, atol=1e-6
    )


def test_BoyerLindquistToCartesianDifferential(bl_differential, cartesian_differential):
    M = 1e24
    a = 0.75
    to_cartesian_differential = bl_differential.cartesian_differential(M=M, a=a)
    assert_allclose(
        to_cartesian_differential.values(),
        cartesian_differential.values(),
        rtol=0.0,
        atol=1e-6,
    )


def test_BoyerLindquistToSphericalDifferential(bl_differential, spherical_differential):
    M = 1e24
    a = 0.75
    to_spherical_differential = bl_differential.spherical_differential(M=M, a=a)
    assert_allclose(
        to_spherical_differential.values(),
        spherical_differential.values(),
        rtol=0.0,
        atol=1e-6,
    )


def test_cycle_BLSphericalDifferential(bl_differential2):
    M = 1e24
    a = 0.75
    bl_diff = bl_differential2
    sph_diff = bl_diff.spherical_differential(M=M, a=a)
    bl_diff2 = sph_diff.bl_differential(M=M, a=a)
    assert_allclose(bl_diff2.values(), bl_diff.values(), rtol=0.0, atol=1e-6)


def test_cycle_BLCartesianDifferential(bl_differential2):
    M = 1e24
    a = 0.75
    bl_diff = bl_differential2
    cart_diff = bl_diff.cartesian_differential(M=M, a=a)
    bl_diff2 = cart_diff.bl_differential(M=M, a=a)
    assert_allclose(bl_diff2.values(), bl_diff.values(), rtol=0.0, atol=1e-6)


def test_state(cartesian_differential):
    M = 1e24
    ms = Schwarzschild(M=M)

    sph = cartesian_differential.spherical_differential()

    st_cart = cartesian_differential.state(metric=ms)
    st_sph = sph.state(metric=ms)

    assert st_cart[0] == st_sph[0]


# Tests for object.__str__ and object.__repr__
@mock.patch("sys.stdout", new_callable=StringIO)
def test_str_diff_objects(
    mock_stdout, cartesian_differential, spherical_differential, bl_differential
):
    print(str(cartesian_differential))
    assert "object at 0x" not in mock_stdout.getvalue()

    print(str(spherical_differential))
    assert "object at 0x" not in mock_stdout.getvalue()

    print(str(bl_differential))
    assert "object at 0x" not in mock_stdout.getvalue()


@mock.patch("sys.stdout", new_callable=StringIO)
def test_repr_diff_objects(
    mock_stdout, cartesian_differential, spherical_differential, bl_differential
):
    print(repr(cartesian_differential))
    assert "object at 0x" not in mock_stdout.getvalue()

    print(repr(spherical_differential))
    assert "object at 0x" not in mock_stdout.getvalue()

    print(repr(bl_differential))
    assert "object at 0x" not in mock_stdout.getvalue()


# Tests for object.velocity()
def test_velocity(cartesian_differential, spherical_differential, bl_differential):
    M = 1e24
    a = 0.
    ms = Schwarzschild(M=M)
    mk = Kerr(coords="BL", M=M, a=a)

    def with_numpy_array(differential_obj, metric):
        assert_allclose(
            differential_obj.values()[4:],
            differential_obj.velocity(metric=metric)[1:],
            1e-10,
            1e-15,
        )

    with_numpy_array(cartesian_differential, ms)
    with_numpy_array(spherical_differential, ms)
    with_numpy_array(bl_differential, mk)


def test_velocity2(cartesian_differential, spherical_differential, bl_differential):
    M = 1e24
    a = 0.
    ms = Schwarzschild(M=M)
    mk = Kerr(coords="BL", M=M, a=a)

    cd, sd, bd = cartesian_differential, spherical_differential, bl_differential
    assert_allclose(cd.velocity(metric=ms)[1:], [cd.v_x, cd.v_y, cd.v_z])
    assert_allclose(sd.velocity(metric=ms)[1:], [sd.v_r, sd.v_th, sd.v_p])
    assert_allclose(bd.velocity(metric=mk)[1:], [bd.v_r, bd.v_th, bd.v_p])
