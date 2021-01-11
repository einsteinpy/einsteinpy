import pytest

from einsteinpy.utils import CoordinateError


def test_CoordinateError_class():
    """
    Tests, if the errors raised with the CoordinateError class \
    ``einsteinpy.coordinates.utils.CoordinateError`` are \
    appropriate

    """
    err_string = " ".join(
        str(i)
        for i in [
            "Error",
            404,
            chr(0xA),
            "\bPage Found =",
            False,
            "\b, Please contact :",
            11.223344,
            None,
            0x1F76,
            12e-4,
        ]
    )
    err = CoordinateError(
        "Error",
        404,
        chr(0xA),
        "\bPage Found =",
        False,
        "\b, Please contact :",
        11.223344,
        None,
        0x1F76,
        12e-4,
    )

    assert hasattr(err, "message")
    assert str(err) == repr(err)
    assert err_string == err.message
