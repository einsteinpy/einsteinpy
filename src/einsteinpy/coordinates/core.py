import numpy as np
from astropy import units as u

from einsteinpy import constant
from einsteinpy.coordinates.conversion import (
    BaseCoordinateConversion,
    BoyerLindquistConversion,
    CartesianConversion,
    SphericalConversion,
)

_c = constant.c.value


class BaseCoordinate(BaseCoordinateConversion):
    """
    Base Coordinate class

    """

    def __init__(self, e0, e1, e2, e3):
        """
        Constructor

        Parameters
        ----------
        e0 : float
            Time
        e1 : float
            x-Component of 3-Position
        e2 : float
            y-Component of 3-Position
        e3 : float
            z-Component of 3-Position

        """
        super().__init__(e0.si.value, e1.si.value, e2.si.value, e3.si.value)
        self.e0 = e0
        self.e1 = e1
        self.e2 = e2
        self.e3 = e3


class Cartesian(BaseCoordinate, CartesianConversion):
    """
    Class for defining 3-Position & 4-Position in Cartesian Coordinates
    using SI units

    """

    @u.quantity_input(e0=u.s, e1=u.m, e2=u.m, e3=u.m)
    def __init__(self, e0, e1, e2, e3):
        super().__init__(e0, e1, e2, e3)
        self.system = "Cartesian"
        self._dimension = {
            "e0": self.e0,
            "e1": self.e1,
            "e2": self.e2,
            "e3": self.e3,
            "system": self.system,
        }
        self._dimension_order = ("e0", "e1", "e2", "e3")

    def __str__(self):
        return f"Cartesian Coordinates: \n \
            e0 = ({self.e0}), e1 = ({self.e1}), e2 = ({self.e2}), e3 = ({self.e3})"

    __repr__ = __str__

    def __getitem__(self, item):
        """
        Method to return coordinates
        Objects are subscriptable with both explicit names of \
        parameters and integer indices

        Parameters
        ----------
        item : str or int
            Name of the parameter or its index
            If ``system`` is provided, while initializing, \
            name of the coordinate is returned

        """
        if isinstance(item, (int, np.integer)):
            return self._dimension[self._dimension_order[item]]
        return self._dimension[item]

    def position(self):
        """
        Returns Position 4-Vector in SI units

        Returns
        -------
        tuple
            4-Tuple, containing Position 4-Vector in SI units

        """
        return (
            _c * self.e0.si.value,
            self.e1.si.value,
            self.e2.si.value,
            self.e3.si.value,
        )

    def to_spherical(self, **kwargs):
        """
        Method for conversion to Spherical Polar Coordinates

        Other Parameters
        ----------------
        **kwargs : dict
            Keyword Arguments

        Returns
        -------
        ~einsteinpy.coordinates.core.Spherical
            Spherical representation of the Cartesian Coordinates

        """
        e0, e1, e2, e3 = self.convert_spherical()

        return Spherical(e0 * u.s, e1 * u.m, e2 * u.rad, e3 * u.rad)

    def to_bl(self, **kwargs):
        """
        Method for conversion to Boyer-Lindquist (BL) Coordinates

        Parameters
        ----------
        **kwargs : dict
            Keyword Arguments
            Expects two arguments, ``M and ``a``, as described below

        Other Parameters
        ----------------
        M : float
            Mass of gravitating body
            Required to calculate ``alpha``, the rotational length \
            parameter
        a : float
            Spin Parameter of gravitating body
            0 <= a <= 1
            Required to calculate ``alpha``, the rotational length \
            parameter

        Returns
        -------
        ~einsteinpy.coordinates.core.BoyerLindquist
            Boyer-Lindquist representation of the Cartesian Coordinates

        """
        M, a = kwargs["M"], kwargs["a"]
        e0, e1, e2, e3 = self.convert_bl(M=M, a=a)

        return BoyerLindquist(e0 * u.s, e1 * u.m, e2 * u.rad, e3 * u.rad)


class Spherical(BaseCoordinate, SphericalConversion):
    """
    Class for defining 3-Position & 4-Position in Spherical Polar Coordinates \
    using SI units

    """

    @u.quantity_input(e0=u.s, e1=u.m, e2=u.rad, e3=u.rad)
    def __init__(self, e0, e1, e2, e3):
        super().__init__(e0, e1, e2, e3)
        self.system = "Spherical"
        self._dimension = {
            "e0": self.e0,
            "e1": self.e1,
            "e2": self.e2,
            "e3": self.e3,
            "system": self.system,
        }
        self._dimension_order = ("e0", "e1", "e2", "e3")

    def __str__(self):
        return f"Spherical Polar Coordinates: \n \
            e0 = ({self.e0}), e1 = ({self.e1}), e2 = ({self.e2}), e3 = ({self.e3})"

    __repr__ = __str__

    def __getitem__(self, item):
        """
        Method to return coordinates
        Objects are subscriptable with both explicit names of \
        parameters and integer indices

        Parameters
        ----------
        item : str or int
            Name of the parameter or its index
            If ``system`` is provided, while initializing, \
            name of the coordinate is returned

        """
        if isinstance(item, (int, np.integer)):
            return self._dimension[self._dimension_order[item]]
        return self._dimension[item]

    def position(self):
        """
        Returns Position 4-Vector in SI units

        Returns
        -------
        tuple :
            4-Tuple, containing Position 4-Vector in SI units

        """
        return (
            _c * self.e0.si.value,
            self.e1.si.value,
            self.e2.si.value,
            self.e3.si.value,
        )

    def to_cartesian(self, **kwargs):
        """
        Method for conversion to Cartesian Coordinates

        Other Parameters
        ----------------
        **kwargs : dict
            Keyword Arguments

        Returns
        -------
        ~einsteinpy.coordinates.core.Cartesian
            Cartesian representation of the Spherical Polar Coordinates

        """
        e0, e1, e2, e3 = self.convert_cartesian()

        return Cartesian(e0 * u.s, e1 * u.m, e2 * u.m, e3 * u.m)

    def to_bl(self, **kwargs):
        """
        Method for conversion to Boyer-Lindquist (BL) Coordinates

        Parameters
        ----------
        **kwargs : dict
            Keyword Arguments
            Expects two arguments, ``M and ``a``, as described below

        Other Parameters
        ----------------
        M : float
            Mass of gravitating body
            Required to calculate ``alpha``, the rotational length \
            parameter
        a : float
            Spin Parameter of gravitating body
            0 <= a <= 1
            Required to calculate ``alpha``, the rotational length \
            parameter

        Returns
        -------
        ~einsteinpy.coordinates.core.BoyerLindquist
            Boyer-Lindquist representation of the Spherical \
            Polar Coordinates

        """
        M, a = kwargs["M"], kwargs["a"]
        e0, e1, e2, e3 = self.convert_bl(M=M, a=a)

        return BoyerLindquist(e0 * u.s, e1 * u.m, e2 * u.rad, e3 * u.rad)


class BoyerLindquist(BaseCoordinate, BoyerLindquistConversion):
    """
    Class for defining 3-Position & 4-Position in Boyer-Lindquist Coordinates \
    using SI units

    """

    @u.quantity_input(e0=u.s, e1=u.m, e2=u.rad, e3=u.rad)
    def __init__(self, e0, e1, e2, e3):
        super().__init__(e0, e1, e2, e3)
        self.system = "BoyerLindquist"
        self._dimension = {
            "e0": self.e0,
            "e1": self.e1,
            "e2": self.e2,
            "e3": self.e3,
            "system": self.system,
        }
        self._dimension_order = ("e0", "e1", "e2", "e3")

    def __str__(self):
        return f"Boyer-Lindquist Coordinates: \n \
            e0 = ({self.e0}), e1 = ({self.e1}), e2 = ({self.e2}), e3 = ({self.e3})"

    __repr__ = __str__

    def __getitem__(self, item):
        """
        Method to return coordinates
        Objects are subscriptable with both explicit names of \
        parameters and integer indices

        Parameters
        ----------
        item : str or int
            Name of the parameter or its index
            If ``system`` is provided, while initializing, \
            name of the coordinate is returned

        """
        if isinstance(item, (int, np.integer)):
            return self._dimension[self._dimension_order[item]]
        return self._dimension[item]

    def position(self):
        """
        Returns Position 4-Vector in SI units

        Returns
        -------
        tuple :
            4-Tuple, containing Position 4-Vector in SI units

        """
        return (
            _c * self.e0.si.value,
            self.e1.si.value,
            self.e2.si.value,
            self.e3.si.value,
        )

    def to_cartesian(self, **kwargs):
        """
        Method for conversion to Cartesian Coordinates

        Parameters
        ----------
        **kwargs : dict
            Keyword Arguments
            Expects two arguments, ``M and ``a``, as described below

        Other Parameters
        ----------------
        M : float
            Mass of gravitating body
            Required to calculate ``alpha``, the rotational length \
            parameter
        a : float
            Spin Parameter of gravitating body
            0 <= a <= 1
            Required to calculate ``alpha``, the rotational length \
            parameter

        Returns
        -------
        ~einsteinpy.coordinates.core.Cartesian
            Cartesian representation of the Boyer-Lindquist Coordinates

        """
        M, a = kwargs["M"], kwargs["a"]
        e0, e1, e2, e3 = self.convert_cartesian(M=M, a=a)

        return Cartesian(e0 * u.s, e1 * u.m, e2 * u.m, e3 * u.m)

    def to_spherical(self, **kwargs):
        """
        Method for conversion to Spherical Polar Coordinates

        Parameters
        ----------
        **kwargs : dict
            Keyword Arguments
            Expects two arguments, ``M and ``a``, as described below

        Other Parameters
        ----------------
        M : float
            Mass of gravitating body
            Required to calculate ``alpha``, the rotational length \
            parameter
        a : float
            Spin Parameter of gravitating body
            0 <= a <= 1
            Required to calculate ``alpha``, the rotational length \
            parameter

        Returns
        -------
        ~einsteinpy.coordinates.core.Spherical
            Spherical Polar representation of the \
            Boyer-Lindquist Coordinates

        """
        M, a = kwargs["M"], kwargs["a"]
        e0, e1, e2, e3 = self.convert_spherical(M=M, a=a)

        return Spherical(e0 * u.s, e1 * u.m, e2 * u.rad, e3 * u.rad)
