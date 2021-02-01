import numpy as np
import pytest
from numpy.testing import assert_allclose

from einsteinpy.geodesic.utils import _P, _sch, _kerr, _kerrnewman


@pytest.mark.parametrize(
    "g, g_prms, q, p, time_like, expected",
    [
        (
            _sch,
            (),
            [0., 2.5, np.pi / 6, np.pi / 2],
            [0.1, 2., 2.],
            True,
            [-0.9167333309092672, 0.1, 2., 2.]
        ),
        (
            _kerr,
            (0.9,),
            [0., 25, np.pi / 2, 0.],
            [0., 0, 2.427],
            False,
            [-0.09333038545829885, 0., 0, 2.427]
        ),
        (
            _kerrnewman,
            (0.5, 0.1,),
            [0., 5.5, np.pi / 4, 0.],
            [0.1, -0.2, -4.],
            True,
            [-1.1220908695019767, 0.1, -0.2, -4.]
        ),
    ],
)
def test_P(g, g_prms, q, p, time_like, expected):
    P = _P(g, g_prms, q, p, time_like)

    assert_allclose(P, expected, atol=1e-8, rtol=1e-8)


@pytest.mark.parametrize(
    "x",
    [
        (
            [0., 2.5, np.pi / 6, np.pi / 2],
        ),
        (
            [0., 25, np.pi / 2, 0.],
        ),
        (
            [0., 5.5, np.pi / 4, 0.],
        ),
    ],
)
def test_metrics(x):
    x = x[0]

    # a = Q = 0.
    s = _sch(x).astype(float)
    k = _kerr(x, 0.).astype(float)
    kn = _kerrnewman(x, 0., 0.).astype(float)
    assert_allclose(s, k, atol=1e-8, rtol=1e-8)
    assert_allclose(k, kn, atol=1e-8, rtol=1e-8)
    assert_allclose(kn, s, atol=1e-8, rtol=1e-8)

    # Non-zero Spin
    k = _kerr(x, 0.4).astype(float)
    kn = _kerrnewman(x, 0.4, 0.).astype(float)
    assert_allclose(k, kn, atol=1e-8, rtol=1e-8)
