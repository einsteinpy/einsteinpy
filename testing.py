"""Testing utilities.
"""
import os.path

import pytest


def test():
    """Initiate einsteinpy testing
    """
    pytest.main([os.path.dirname(os.path.abspath(__file__))])
