import numpy as np
from numpy.testing import assert_allclose
from unittest.mock import Mock

from einsteinpy.examples import precession


partial_lambdas = np.array([0.0000e+00, 5.0000e-01, 1.0000e+00, 1.9990e+03, 1.9995e+03, 2.0000e+03])
partial_vecs = np.array(
    [[ 4.00000000e+01,  0.00000000e+00,  2.44929360e-15, 0.00000000e+00,  0.00000000e+00,  3.83405000e+00],
    [ 3.99999197e+01, -4.79255517e-02,  2.44929044e-15, 2.17126283e-04, -2.81284902e-19,  3.83405000e+00],
    [ 3.99996789e+01, -9.58508493e-02,  2.44928097e-15, 4.34253495e-04, -5.62571255e-19,  3.83405000e+00],
    [-6.00084085e+00, -2.52448871e+00,  3.98636802e-16, -2.33673001e-01, -9.49738949e-15,  3.83405000e+00],
    [-6.18465627e+00, -2.28046143e+00,  4.03625070e-16, -2.36036036e-01, -9.50787744e-15,  3.83405000e+00],
    [-6.35745272e+00, -2.03236148e+00,  4.08689510e-16, -2.38256302e-01, -9.51810736e-15,  3.83405000e+00]]
)

precession = Mock()
precession.return_value = (partial_lambdas, partial_vecs)

def test_precession_attr():
    """
    Checks various attributes of precession()

    """
    p = precession()
    L = p[1][0, 5]

    assert p[0].shape[0] == p[1].shape[1]
    assert p[1].shape[1] == 6
    assert_allclose(p[1][:, 4], 0, atol=1e-12, rtol=1e-12)
    assert_allclose(p[1][:, 5], L, atol=1e-12, rtol=1e-12)
