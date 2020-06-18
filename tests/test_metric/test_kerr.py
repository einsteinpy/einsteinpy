import collections

import numpy as np
import pytest
from numpy.testing import assert_allclose

from einsteinpy.metric import Kerr


def test_nonzero_christoffels():
    """
    Compares count of non-zero Christoffel Symbols List in BL coordinates \
    with that generated algorithmically 

    """
    l1 = Kerr.nonzero_christoffels()
    l2 = Kerr.nonzero_christoffels_list_bl
    
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
        l = t % 4
        for m in range(4):
            chl2[i, k, l] += g_contra[i, m] * (
                dgdx[l, m, k] + dgdx[k, m, l] - dgdx[m, k, l]
            )
    chl2 = np.multiply(chl2, 0.5)
    
    assert_allclose(chl2, chl1, rtol=1e-8)


def test_kerr_private_attr():
    """
    Tries to invoke private methods of Kerr

    em_potential_covariant = _private
    em_potential_contravariant = _private
    em_tensor_covariant = _private
    em_tensor_contravariant = _private
    
    """
    obj = Kerr(coords="BL", M=6e24, a=0.8)

    try:
        obj.em_tensor_covariant()
        obj.em_tensor_contravariant()
        obj.em_potential_covariant()
        obj.em_potential_contravariant()
        assert False

    except AttributeError:
        assert True
