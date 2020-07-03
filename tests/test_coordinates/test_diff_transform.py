import sys
from io import StringIO
from unittest import mock

import astropy.units as u
import numpy as np
import pytest
from numpy.testing import assert_allclose

from einsteinpy.coordinates import (
    BoyerLindquistDifferential,
    CartesianDifferential,
    SphericalDifferential
)
from einsteinpy.metric import Schwarzschild, Kerr
from einsteinpy import constant

_c = constant.c.value


@pytest.fixture()
def cartesian_differential():
    return CartesianDifferential(
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
    return SphericalDifferential(
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
    return BoyerLindquistDifferential(
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
    return BoyerLindquistDifferential(
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
    M = 1e24 * u.kg
    a = 0.75 * u.one
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
    M = 1e24 * u.kg
    a = 0.75 * u.one
    to_bl_differential = spherical_differential.bl_differential(M=M, a=a)
    assert_allclose(
        to_bl_differential.values(), bl_differential.values(), rtol=0.0, atol=1e-6
    )


def test_BoyerLindquistToCartesianDifferential(bl_differential, cartesian_differential):
    M = 1e24 * u.kg
    a = 0.75 * u.one
    to_cartesian_differential = bl_differential.cartesian_differential(M=M, a=a)
    assert_allclose(
        to_cartesian_differential.values(),
        cartesian_differential.values(),
        rtol=0.0,
        atol=1e-6,
    )


def test_BoyerLindquistToSphericalDifferential(bl_differential, spherical_differential):
    M = 1e24 * u.kg
    a = 0.75 * u.one
    to_spherical_differential = bl_differential.spherical_differential(M=M, a=a)
    assert_allclose(
        to_spherical_differential.values(),
        spherical_differential.values(),
        rtol=0.0,
        atol=1e-6,
    )


def test_cycle_BLSphericalDifferential(bl_differential2):
    M = 1e24 * u.kg
    a = 0.75 * u.one
    bl_diff = bl_differential2
    sph_diff = bl_diff.spherical_differential(M=M, a=a)
    bl_diff2 = sph_diff.bl_differential(M=M, a=a)
    assert_allclose(bl_diff2.values(), bl_diff.values(), rtol=0.0, atol=1e-6)


def test_cycle_BLCartesianDifferential(bl_differential2):
    M = 1e24 * u.kg
    a = 0.75 * u.one
    bl_diff = bl_differential2
    cart_diff = bl_diff.cartesian_differential(M=M, a=a)
    bl_diff2 = cart_diff.bl_differential(M=M, a=a)
    assert_allclose(bl_diff2.values(), bl_diff.values(), rtol=0.0, atol=1e-6)


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


def test_position(cartesian_differential, spherical_differential, bl_differential):
    cart_pos = [_c * 1e3, 10e3 / np.sqrt(2), 10e3 / np.sqrt(2), 0]
    sph_pos = [_c * 1e3, 10e3, 1.5707963267948966, 0.7853981633974483]
    bl_pos = [_c * 1e3, 10e3, 1.5707963267948966, 0.7853981633974483]

    assert_allclose(cartesian_differential.position(), cart_pos)
    assert_allclose(spherical_differential.position(), sph_pos)
    assert_allclose(bl_differential.position(), bl_pos)


# Tests for object.velocity()
def test_velocity(spherical_differential, bl_differential):
    M = 1e24 * u.kg
    a = 0. * u.one
    ms = Schwarzschild(coords=spherical_differential, M=M)
    mk = Kerr(coords=bl_differential, M=M, a=a)

    def with_numpy_array(differential_obj, metric):
        assert_allclose(
            differential_obj.values()[4:],
            differential_obj.velocity(metric=metric)[1:],
            1e-10,
            1e-15,
        )

    with_numpy_array(spherical_differential, ms)
    with_numpy_array(bl_differential, mk)


def test_velocity2(spherical_differential, bl_differential):
    M = 1e24 * u.kg
    a = 0. * u.one
    ms = Schwarzschild(coords=spherical_differential, M=M)
    mk = Kerr(coords=bl_differential, M=M, a=a)

    sd, bd = spherical_differential, bl_differential
    assert_allclose(sd.velocity(metric=ms)[1:], [sd.v_r.value, sd.v_th.value, sd.v_p.value])
    assert_allclose(bd.velocity(metric=mk)[1:], [bd.v_r.value, bd.v_th.value, bd.v_p.value])


def test_v_t_raises_TypeError(cartesian_differential, spherical_differential, bl_differential):
    M = 1e24 * u.kg
    a = 0. * u.one
    ms = Schwarzschild(coords=spherical_differential, M=M)
    mk = Kerr(coords=bl_differential, M=M, a=a)

    cd, sd, bd = cartesian_differential, spherical_differential, bl_differential

    def cd_s(cd, ms):
        return cd.velocity(metric=ms)

    def bd_s(bd, ms):
        return bd.velocity(metric=ms)

    def cd_k(cd, mk):
        return cd.velocity(metric=mk)

    def sd_k(sd, mk):
        return sd.velocity(metric=mk)

    pytest.raises(TypeError, cd_s, (cd, ms))
    pytest.raises(TypeError, bd_s, (bd, ms))
    pytest.raises(TypeError, cd_k, (cd, mk))
    pytest.raises(TypeError, sd_k, (sd, mk))


def test_v_t_getter_setter(cartesian_differential, spherical_differential, bl_differential):
    M = 1e24 * u.kg
    a = 0. * u.one
    ms = Schwarzschild(coords=spherical_differential, M=M)
    mk = Kerr(coords=bl_differential, M=M, a=a)

    def cd_vt(cartesian_differential, ms):
        cartesian_differential.v_t = (ms,)

    def sd_vt(spherical_differential, mk):
        spherical_differential.v_t = (mk,)

    def bd_vt(bl_differential, ms):
        bl_differential.v_t = (ms,)

    # These 2 should not raise TypeError
    try:
        spherical_differential.v_t = (ms,)
    except TypeError:
        pytest.fail("Unexpected TypeError!")

    try:
        bl_differential.v_t = (mk,)
    except TypeError:
        pytest.fail("Unexpected TypeError!")

    # These 3 should raise TypeError
    pytest.raises(TypeError, cd_vt, (cartesian_differential, ms))
    pytest.raises(TypeError, sd_vt, (spherical_differential, mk))
    pytest.raises(TypeError, bd_vt, (bl_differential, ms))
