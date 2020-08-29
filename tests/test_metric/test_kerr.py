import collections

import pytest
import astropy.units as u
import numpy as np
from numpy.testing import assert_allclose

from einsteinpy.coordinates import BoyerLindquistDifferential
from einsteinpy.metric import Kerr


def test_nonzero_christoffels():
    """
    Compares count of non-zero Christoffel Symbols List in BL coordinates \
    with that generated algorithmically

    """
    mk = Kerr(coords="BL", M=1. * u.kg, a=0.5 * u.one)
    l1 = mk.nonzero_christoffels()
    l2 = mk._nonzero_christoffels_list_bl

    assert collections.Counter(l1) == collections.Counter(l2)


@pytest.fixture
def bl():
    return BoyerLindquistDifferential(
        t=0. * u.s,
        r=0.1 * u.m,
        theta=4 * np.pi / 5 * u.rad,
        phi=0. * u.rad,
        v_r=0. * u.m / u.s,
        v_th=0. * u.rad / u.s,
        v_p=0. * u.rad / u.s
    )


def test_christoffels(bl):
    """
    Compares output produced by optimized function, with that, produced via general method (formula)

    """
    r, theta = 100.0 * u.m, np.pi / 5 * u.rad
    M, a = 6.73317655e26 * u.kg, 0.2 * u.one

    x_vec = bl.position()

    # Output produced by the optimized function
    mk = Kerr(coords=bl, M=M, a=a)
    chl1 = mk.christoffels(x_vec)

    # Calculated using formula
    g_contra = mk.metric_contravariant(x_vec)
    dgdx = mk._dg_dx_bl(x_vec)
    chl2 = np.zeros(shape=(4, 4, 4), dtype=float)
    tmp = np.array([i for i in range(4 ** 3)])
    for t in tmp:
        i = int(t / (4 ** 2)) % 4
        k = int(t / 4) % 4
        index = t % 4
        for m in range(4):
            chl2[i, k, index] += g_contra[i, m] * (
                dgdx[index, m, k] + dgdx[k, m, index] - dgdx[m, k, index]
            )
    chl2 = np.multiply(chl2, 0.5)

    assert_allclose(chl2, chl1, rtol=1e-8)


def test_f_vec_bl_kerr():
    M, a = 6.73317655e26 * u.kg, 0.2 * u.one
    bl = BoyerLindquistDifferential(
        t=0. * u.s,
        r=1e6 * u.m,
        theta=4 * np.pi / 5 * u.rad,
        phi=0. * u.rad,
        v_r=0. * u.m / u.s,
        v_th=0. * u.rad / u.s,
        v_p=2e6 * u.rad / u.s
    )
    f_vec_expected = np.array(
        [3.92128321e+03,  0.00000000e+00,  0.00000000e+00,  2.00000000e+06, -0.00000000e+00,  1.38196394e+18, -1.90211303e+12, -0.00000000e+00]
    )

    mk = Kerr(coords=bl, M=M, a=a)
    state = np.hstack((bl.position(), bl.velocity(mk)))

    f_vec = mk._f_vec(0., state)

    assert isinstance(f_vec, np.ndarray)
    assert_allclose(f_vec_expected, f_vec, rtol=1e-8)
