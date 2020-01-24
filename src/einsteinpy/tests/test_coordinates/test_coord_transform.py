import sys
from io import StringIO

import astropy.units as u
import numpy as np
import pytest
from numpy.testing import assert_allclose

from einsteinpy.coordinates import BoyerLindquist, Cartesian, Spherical


@pytest.fixture()
def cartesian():
    return Cartesian(20.0 * u.km, 311e3 * u.m, 210e5 * u.cm)


@pytest.fixture()
def spherical():
    return Spherical(
        375.79382645275 * u.km, 0.9778376650369 * u.rad, 1.5065760775947 * u.rad
    )


@pytest.fixture()
def boyerlindquist():
    return BoyerLindquist(
        375.79382645275 * u.km,
        0.9778376650369 * u.rad,
        1.5065760775947 * u.rad,
        0.0 * u.cm,
    )


def test_CartesianToSpherical(cartesian, spherical):
    to_spherical = cartesian.to_spherical()
    assert_allclose(
        to_spherical.si_values(), spherical.si_values(), rtol=0.0, atol=1e-6
    )


def test_CartesianToBoyerLindquist(cartesian, boyerlindquist):
    bl = boyerlindquist
    to_bl = cartesian.to_bl(bl.a)
    assert_allclose(to_bl.si_values(), bl.si_values(), rtol=0.0, atol=1e-6)


def test_SphericalToCartesian(spherical, cartesian):
    to_cartesian = spherical.to_cartesian()
    assert_allclose(
        to_cartesian.si_values(), cartesian.si_values(), rtol=0.0, atol=1e-6
    )


def test_SphericalToBoyerLindquist(spherical, boyerlindquist):
    bl = boyerlindquist
    to_bl = spherical.to_bl(bl.a)
    assert_allclose(to_bl.si_values(), bl.si_values(), rtol=0.0, atol=1e-6)


def test_BoyerLindquistToCartesian(boyerlindquist, cartesian):
    bl = boyerlindquist
    to_cartesian = bl.to_cartesian()
    assert_allclose(
        to_cartesian.si_values(), cartesian.si_values(), rtol=0.0, atol=1e-6
    )


def test_BoyerLindquistToSpherical(boyerlindquist, spherical):
    bl = boyerlindquist
    to_spherical = bl.to_spherical()
    assert_allclose(
        to_spherical.si_values(), spherical.si_values(), rtol=0.0, atol=1e-6
    )


@pytest.mark.parametrize(
    "cart, a",
    [
        (Cartesian(10 * u.m, 10 * u.m, 0 * u.m), 0.7 * u.m),
        (Cartesian(-732.0 * u.km, 456 * u.km, -90 * u.km), 9.0 * u.km),
        (Cartesian(-0.732 * u.km, -1.456 * u.km, 90 * u.m), 21 * u.m),
        (Cartesian(2 * u.m, -1 * u.m, 0 * u.cm), 0 * u.cm),
    ],
)
def test_cycle_CartesianBL(cart, a):
    bl = cart.to_bl(a)
    cart2 = bl.to_cartesian()
    assert_allclose(cart2.si_values(), cart.si_values(), rtol=0.0, atol=1e-6)


@pytest.mark.parametrize(
    "bl",
    [
        BoyerLindquist(3 * u.m, 120 * u.deg, 175 * u.deg, 0.3 * u.m),
        BoyerLindquist(100 * u.m, 1e-3 * u.rad, 0 * u.deg, 6.969 * u.m),
    ],
)
def test_cycle_BLSpherical(bl):
    sph = bl.to_spherical()
    bl2 = sph.to_bl(bl.a)
    assert_allclose(bl2.si_values(), bl.si_values(), rtol=0.0, atol=1e-6)


def test_cartesian_norm():
    test_data_x = 1 * u.km
    test_data_y = 1 * u.km
    test_data_z = 1 * u.km

    test_obj = Cartesian(x=test_data_x, y=test_data_y, z=test_data_z)
    assert_allclose(
        (test_obj.norm()).si.value, (np.sqrt(3) * u.km).si.value, rtol=0, atol=1e-5
    )


def test_cartesian_dot():
    test_data_x = 3 * u.km
    test_data_y = 3 * u.km
    test_data_z = 3 * u.km

    test_data_x_target = 3 * u.km
    test_data_y_target = 3 * u.km
    test_data_z_target = 3 * u.km

    test_obj = Cartesian(x=test_data_x, y=test_data_y, z=test_data_z)
    test_target_obj = Cartesian(
        x=test_data_x_target, y=test_data_y_target, z=test_data_z_target
    )

    assert_allclose(
        (test_obj.dot(test_target_obj)).si.value,
        (27 * u.km * u.km).si.value,
        rtol=0,
        atol=1e-5,
    )


# Tests for object.__repr__ and object.__str__


def test_print_core_objects(cartesian, spherical, boyerlindquist):
    old_stdout = sys.stdout
    # change stdout
    result = StringIO()
    sys.stdout = result

    print(str(cartesian))
    assert "object at 0x" not in result.getvalue()

    print(str(spherical))
    assert "object at 0x" not in result.getvalue()

    print(str(boyerlindquist))
    assert "object at 0x" not in result.getvalue()

    # again switch to old stdout
    sys.stdout = old_stdout
