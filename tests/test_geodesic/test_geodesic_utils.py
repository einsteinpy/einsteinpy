import warnings

import numpy as np
import pytest
from numpy.testing import assert_allclose
from scipy.optimize import fsolve

from einsteinpy.geodesic.utils import _energy, _sphToCart


@pytest.mark.parametrize(
    "q, p, a, mu, expected",
    [
        (
            [2.5, np.pi / 6, np.pi / 2],
            [0.1, 2., 2.],
            0,
            -1,
            0.9167333309092672
        ),
        (
            [25, np.pi / 2, 0.],
            [0., 0, 2.427],
            0.99,
            -1,
            0.9639802199360611
        ),
        (
            [5.5, np.pi / 4, 0.],
            [0.1, -0.2, -4.],
            0.99,
            0.,
            0.7681630092559442
        ),
    ],
)
def test_energy(q, p, a, mu, expected):
    E = fsolve(_energy, 1.0, args=(q, p, a, mu))[0]

    assert_allclose(E, expected, atol=1e-8, rtol=1e-8)


@pytest.mark.parametrize(
    "q, expected",
    [
        (
            [2.5, np.pi / 6, np.pi / 2],
            [7.654042494670957e-17, 1.2499999999999998, 2.165063509461097]
        ),
        (
            [25, np.pi / 2, 0.],
            [25.0, 0.0, 1.5308084989341915e-15]
        ),
        (
            [5.5, np.pi / 4, 0.],
            [3.8890872965260117, 0.0, 3.8890872965260117]
        ),
    ],
)
def test_sphToCart(q, expected):
    qc = _sphToCart(*q)

    assert_allclose(qc, expected, atol=1e-8, rtol=1e-8)
