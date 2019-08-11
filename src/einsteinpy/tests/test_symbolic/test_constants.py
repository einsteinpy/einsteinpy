from einsteinpy.symbolic import constants, get_constant


def test_get_constant():
    c1 = constants.c
    c2 = get_constant("c")
    assert c1 - c2 == 0


def test_SymbolicConstant_descriptive_name():
    ld = constants.Cosmo_Const
    assert ld._descriptive_name == ld.descriptive_name
