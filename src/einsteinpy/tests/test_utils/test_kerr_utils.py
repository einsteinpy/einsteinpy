import collections

import pytest

from einsteinpy.utils import kerr_utils


def test_nonzero_christoffels():
    l1 = kerr_utils.nonzero_christoffels()
    l2 = kerr_utils.nonzero_christoffels_list
    assert collections.Counter(l1) == collections.Counter(l2)
