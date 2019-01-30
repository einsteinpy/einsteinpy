


def test_primes():
    from ..example_mod import primes
    assert primes(10) == [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]


def test_deprecation():
    import warnings
    warnings.warn(
        "This is deprecated, but shouldn't raise an exception, unless "
        "enable_deprecations_as_exceptions() called from conftest.py",
        DeprecationWarning)
