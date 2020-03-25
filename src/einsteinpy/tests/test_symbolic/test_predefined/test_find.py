import pytest

from einsteinpy.symbolic.predefined import find


def test_find_correct():
    names = ["AntiDeSitter", "AntiDeSitterStatic", "DeSitter"]
    find_list = find("sit")
    assert set(names).issubset(set(find_list))


def test_find_incorrect():
    assert find("doesnotexist") == []
    assert find("find") == []
