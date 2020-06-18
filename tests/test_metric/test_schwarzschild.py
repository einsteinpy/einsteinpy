import pytest

from einsteinpy.metric import Schwarzschild


def test_schwarzschild_private_attr():
    """
    Tries to invoke private methods of Kerr

    em_potential_covariant = _private
    em_potential_contravariant = _private
    em_tensor_covariant = _private
    em_tensor_contravariant = _private

    """
    obj = Schwarzschild(M=6e24)

    try:
        obj.em_tensor_covariant()
        obj.em_tensor_contravariant()
        obj.em_potential_covariant()
        obj.em_potential_contravariant()
        assert False

    except AttributeError:
        assert True
