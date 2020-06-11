import pytest

from einsteinpy.symbolic.predefined import find


def test_find_correct():
    names = ["AntiDeSitter", "AntiDeSitterStatic", "DeSitter"]
    namesMin = ["Minkowski", "MinkowskiCartesian", "MinkowskiPolar"]
    find_list = find("sit")
    find_listMin = find("min")
    assert set(names).issubset(set(find_list))
    assert set(namesMin).issubset(set(find_listMin))


def test_find_incorrect():
    assert find("doesnotexist") == []
    assert find("find") == []
