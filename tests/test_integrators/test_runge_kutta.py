import warnings

import numpy as np
import pytest

from einsteinpy import integrators


def test_RK45_RuntimeWarning():
    def fun(t, y):
        return 2 * y

    cl = integrators.RK45(fun, 3.0, np.array([3, 2, 1.0]), 8, 2.0)
    with warnings.catch_warnings(record=True) as w:
        for i in range(10):
            cl.step()
        assert len(w) == 1 and issubclass(w[-1].category, RuntimeWarning)
