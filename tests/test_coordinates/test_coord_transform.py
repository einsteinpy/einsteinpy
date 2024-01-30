import sys
from io import StringIO
from unittest import mock

import astropy.units as u
import numpy as np
import pytest
from numpy.testing import assert_allclose

from einsteinpy.coordinates import BoyerLindquist, Cartesian, Spherical
from einsteinpy import constant

_c = constant.c.value


@pytest.fixture()
def cartesian():
    return Cartesian(
        0.0 * u.s,
        2e4 * u.m,
        311e3 * u.m,
        210e3 * u.m
    )


@pytest.fixture()
def spherical():
    return Spherical(
        0.0 * u.s,
        375793.82645275 * u.m,
        0.9778376650369 * u.rad,
        1.5065760775947 * u.rad
    )


@pytest.fixture()
def boyerlindquist():
    return BoyerLindquist(
        0.0 * u.s,
        375793.82645275 * u.m,
        0.9778376650369 * u.rad,
        1.5065760775947 * u.rad
    )


def test_CartesianToSpherical(cartesian, spherical):
    to_spherical = cartesian.to_spherical()
    assert_allclose(
        to_spherical.values(), spherical.values(), rtol=0.0, atol=1e-6
    )


def test_CartesianToBoyerLindquist(cartesian, boyerlindquist):
    M = 1e24 * u.kg
    a = 0.75 * u.one
    bl = boyerlindquist
    to_bl = cartesian.to_bl(M=M, a=a)
    assert_allclose(to_bl.values(), bl.values(), rtol=0.0, atol=1e-6)


def test_SphericalToCartesian(spherical, cartesian):
    to_cartesian = spherical.to_cartesian()
    assert_allclose(
        to_cartesian.values(), cartesian.values(), rtol=0.0, atol=1e-6
    )


def test_SphericalToBoyerLindquist(spherical, boyerlindquist):
    M = 1e24 * u.kg
    a = 0.75 * u.one
    bl = boyerlindquist
    to_bl = spherical.to_bl(M=M, a=a)
    assert_allclose(to_bl.values(), bl.values(), rtol=0.0, atol=1e-6)


def test_BoyerLindquistToCartesian(boyerlindquist, cartesian):
    M = 1e24 * u.kg
    a = 0.75 * u.one
    bl = boyerlindquist
    to_cartesian = bl.to_cartesian(M=M, a=a)
    assert_allclose(
        to_cartesian.values(), cartesian.values(), rtol=0.0, atol=1e-6
    )


def test_BoyerLindquistToSpherical(boyerlindquist, spherical):
    M = 1e24 * u.kg
    a = 0.75 * u.one
    bl = boyerlindquist
    to_spherical = bl.to_spherical(M=M, a=a)
    assert_allclose(
        to_spherical.values(), spherical.values(), rtol=0.0, atol=1e-6
    )


@pytest.mark.parametrize(
    "cart, M, a",
    [
        (
            Cartesian(1e0 * u.s, 10 * u.m, 10 * u.m, 0 * u.m),
            1e24 * u.kg,
            0.7 * u.one
        ),
        (
            Cartesian(1e1 * u.s, -732e2 * u.m, 456e2 * u.m, -90e2 * u.m),
            3e10 * u.kg,
            0.9 * u.one
        ),
        (
            Cartesian(1e2 * u.s, -0.732e2 * u.m, -1.456e2 * u.m, 90 * u.m),
            55e20 * u.kg,
            0.3 * u.one
        ),
        (
            Cartesian(1e3 * u.s, 2 * u.m, -1 * u.m, 0 * u.m),
            7e10 * u.kg,
            0 * u.one
        ),
    ],
)
def test_cycle_CartesianBL(cart, M, a):
    bl = cart.to_bl(M=M, a=a)
    cart2 = bl.to_cartesian(M=M, a=a)
    assert_allclose(cart2.values(), cart.values(), rtol=0.0, atol=1e-6)


@pytest.mark.parametrize(
    "bl, M, a",
    [
        (
            BoyerLindquist(1e2 * u.s, 3 * u.m, 2 * np.pi / 3 * u.rad, 3.0543261909900767 * u.rad),
            1e24 * u.kg,
            0.3 * u.one
        ),
        (
            BoyerLindquist(1e3 * u.s, 100 * u.m, 1e-3 * u.rad, 0 * u.rad),
            4e20 * u.kg,
            0.9 * u.one
        )
    ],
)
def test_cycle_BLSpherical(bl, M, a):
    sph = bl.to_spherical(M=M, a=a)
    bl2 = sph.to_bl(M=M, a=a)
    assert_allclose(bl2.values(), bl.values(), rtol=0.0, atol=1e-6)


def test_position(cartesian, spherical, boyerlindquist):
    cart_pos = np.array([_c * cartesian[0].value, cartesian[1].value, cartesian[2].value, cartesian[3].value])
    sph_pos = np.array([_c * spherical[0].value, spherical[1].value, spherical[2].value, spherical[3].value])
    bl_pos = np.array([
        _c * boyerlindquist[0].value,
        boyerlindquist[1].value,
        boyerlindquist[2].value,
        boyerlindquist[3].value
    ])

    assert_allclose(cartesian.position(), cart_pos)
    assert_allclose(spherical.position(), sph_pos)
    assert_allclose(boyerlindquist.position(), bl_pos)


# Tests for object.__str__ and object.__repr__
@mock.patch("sys.stdout", new_callable=StringIO)
def test_str_core_objects(mock_stdout, cartesian, spherical, boyerlindquist):
    print(str(cartesian))
    assert "object at 0x" not in mock_stdout.getvalue()

    print(str(spherical))
    assert "object at 0x" not in mock_stdout.getvalue()

    print(str(boyerlindquist))
    assert "object at 0x" not in mock_stdout.getvalue()


@mock.patch("sys.stdout", new_callable=StringIO)
def test_repr_core_objects(mock_stdout, cartesian, spherical, boyerlindquist):
    print(repr(cartesian))
    assert "object at 0x" not in mock_stdout.getvalue()

    print(repr(spherical))
    assert "object at 0x" not in mock_stdout.getvalue()

    print(repr(boyerlindquist))
    assert "object at 0x" not in mock_stdout.getvalue()


def test_coordinate_subscripting(cartesian, spherical, boyerlindquist):
    assert cartesian["e1"] == cartesian.e1 == cartesian[1] == cartesian[-3]
    assert spherical["e2"] == spherical.e2 == spherical[2] == spherical[-2]
    assert (
        boyerlindquist["e3"] == 
        boyerlindquist.e3 == 
        boyerlindquist[3] == 
        boyerlindquist[-1]
    )
