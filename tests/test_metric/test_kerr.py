import collections

import pytest
import numpy as np
from numpy.testing import assert_allclose

from einsteinpy.metric import Kerr


def test_nonzero_christoffels():
    """
    Compares count of non-zero Christoffel Symbols List in BL coordinates \
    with that generated algorithmically

    """
    mk = Kerr(coords="BL", M=1., a=0.5)
    l1 = mk.nonzero_christoffels()
    l2 = mk._nonzero_christoffels_list_bl

    assert collections.Counter(l1) == collections.Counter(l2)


def test_christoffels():
    """
    Compares output produced by optimized function, with that, produced via general method (formula)

    """
    r, theta = 100.0, np.pi / 5
    M, a = 6.73317655e26, 0.2
    x_vec = np.array([0., r, theta, 0.])

    # Output produced by the optimized function
    mk = Kerr(coords="BL", M=M, a=a)
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


@pytest.mark.parametrize(
    "func_ks",
    [
        (
            Kerr(coords="KS", M=1e22, a=0.7).metric_covariant
        ),
        (
            Kerr(coords="KS", M=1e22, a=0.7).metric_contravariant
        ),
        (
            Kerr(coords="KS", M=1e22, a=0.7).christoffels
        ),
        (
            Kerr(coords="KS", M=1e22, a=0.7)._dg_dx_ks
        ),
    ],
)
def test_ks_raises_NotImplementedError(func_ks):
    """
    Tests, if NotImplementedError is raised, when Kerr-Schild coordinates are used

    """
    x_vec = np.array([0., 5.5, 2 * np.pi / 5, 0.])

    try:
        func_ks(x_vec)
        assert False

    except NotImplementedError:
        assert True


def test_f_vec_ks_raises_NotImplementedError():
    """
    Tests, if NotImplementedError is raised by ``f_vec_ks()``, when Kerr-Schild coordinates are used

    """
    x_vec = np.array([0., 5.5, 2 * np.pi / 5, 0.])

    try:
        Kerr(coords="KS", M=1e22, a=0.7).f_vec(0., x_vec)
        assert False

    except NotImplementedError:
        assert True


def test_singularities_ks_raises_NotImplementedError():
    """
    Tests, if a NotImplementedError is raised, when KerrSchild coordinates \
    are used with ``singularities()``

    """
    mk = Kerr(coords="KS", M=1e22, a=0.5)

    try:
        mk_sing = mk.singularities()
        assert False

    except NotImplementedError:
        assert True
