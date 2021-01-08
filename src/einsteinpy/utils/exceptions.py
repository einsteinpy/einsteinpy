"""Docstring for exceptions.py module

This module defines the ``BaseError`` class which is the base class for all \
custom Errors in EinsteinPy, and the ``CoordinateError`` class, which is a child class \
used for raising exceptions when the geometry does not support \
the supplied coordinate system.

"""


class BaseError(Exception):
    """
    Base class for custom errors

    """

    def __init__(self, *args, **kwargs):
        """
        | Constructor
        | Joins ``args`` into a ``message`` string

        Parameters
        ------------
        *args : iterable
            Other arguments
        **kwargs : dict
            Keyword arguments

        """
        self.message = " ".join(str(i) for i in args)

    def __repr__(self):
        return f"{self.__class__.__name__} : {self.message}"

    def __str__(self):
        return self.__repr__()


class CoordinateError(BaseError):
    """
    Error class for invalid coordinate operations

    """

    def __init__(self, *args, **kwargs):
        """
        | Constructor
        | Joins ``args`` into a ``message`` string

        Parameters
        ------------
        *args : iterable
            Other arguments
        **kwargs : dict
            Keyword arguments

        """
        super().__init__(*args, **kwargs)
