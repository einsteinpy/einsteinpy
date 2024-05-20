import numpy as np
from astropy import units as u

from einsteinpy import constant, metric
from einsteinpy.coordinates.conversion import (
    BaseCoordinateConversion,
    BoyerLindquistConversion,
    CartesianConversion,
    SphericalConversion,
)
from einsteinpy.coordinates.utils import v0
from einsteinpy.utils import CoordinateError

_c = constant.c.value


class BaseDifferential(BaseCoordinateConversion):
    """
    Base Differential Class

    """

    def __init__(self, e0, e1, e2, e3, u0, u1, u2):
        """
        Constructor

        Parameters
        ----------
        e0 : ~astropy.units.quantity.Quantity
            Time
        e1 : ~astropy.units.quantity.Quantity
            x-Component of 3-Position
        e2 : ~astropy.units.quantity.Quantity
            y-Component of 3-Position
        e3 : ~astropy.units.quantity.Quantity
            z-Component of 3-Position
        u0 : ~astropy.units.quantity.Quantity, optional
            x-Component of 3-Velocity
        u1 : ~astropy.units.quantity.Quantity, optional
            y-Component of 3-Velocity
        u2 : ~astropy.units.quantity.Quantity, optional
            z-Component of 3-Velocity

        """
        super().__init__(
            e0.si.value,
            e1.si.value,
            e2.si.value,
            e3.si.value,
            u0.si.value,
            u1.si.value,
            u2.si.value,
        )
        self.e0 = e0
        self.e1 = e1
        self.e2 = e2
        self.e3 = e3
        self._v_t = None
        self.u0 = u0
        self.u1 = u1
        self.u2 = u2

    def __str__(self):
        return f"{self.__class__.__name__} Coordinates: \n\
            e0 = ({self.e0}), e1 = ({self.e1}), e2 = ({self.e2}), e3 = ({self.e3})\n\
            v_t: {self.v_t}, u0: {self.u0}, u1: {self.u1}, u2: {self.u2}"

    __repr__ = __str__

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

    def velocity(self, metric):
        """
        Returns Velocity 4-Vector in SI units

        Parameters
        ----------
        metric : ~einsteinpy.metric.*
            Metric object, in which the coordinates are defined

        Returns
        -------
        tuple
            4-Tuple, containing Velocity 4-Vector in SI units

        """
        # Setting _v_t
        self.v_t = (metric,)

        return (
            self._v_t.value,
            self.u0.si.value,
            self.u1.si.value,
            self.u2.si.value,
        )


class CartesianDifferential(BaseDifferential, CartesianConversion):
    """
    Class for defining 3-Velocity & 4-Velocity in Cartesian Coordinates \
    using SI units

    """

    @u.quantity_input(
        e0=u.s, e1=u.m, e2=u.m, e3=u.m, u0=u.m / u.s, u1=u.m / u.s, u2=u.m / u.s
    )
    def __init__(self, e0, e1, e2, e3, u0, u1, u2):
        super().__init__(e0, e1, e2, e3, u0, u1, u2)
        self.system = "Cartesian"

    @property
    def v_t(self):
        """
        Returns the Timelike component of 4-Velocity

        """
        return self._v_t

    @v_t.setter
    def v_t(self, args):
        """
        Sets the value of the Time-like component of 4-Velocity

        Parameters
        ----------
        args : tuple
            1-tuple containing the ~einsteinpy.metric.* object, \
            in which the coordinates are defined

        Raises
        ------
        CoordinateError
            If ``metric`` object has been instantiated with a coordinate system, \
            other than Cartesian Coordinates.

        """
        g = args[0]
        if self.system != g.coords.system:
            raise CoordinateError(
                f"Metric object has been instantiated with a coordinate system, ( {g.coords.system} )"
                " other than Cartesian Coordinates."
            )

        g_cov_mat = g.metric_covariant(self.position())

        v_t = v0(g_cov_mat, self.u0.si.value, self.u1.si.value, self.u2.si.value)

        self._v_t = v_t * u.m / u.s

    def spherical_differential(self, **kwargs):
        """
        Converts to Spherical Polar Coordinates

        Parameters
        ----------
        **kwargs : dict
            Keyword Arguments

        Returns
        -------
        ~einsteinpy.coordinates.differential.SphericalDifferential
            Spherical Polar representation of velocity

        """
        e0, e1, e2, e3, u0, u1, u2 = self.convert_spherical()
        return SphericalDifferential(
            e0 * u.s,
            e1 * u.m,
            e2 * u.rad,
            e3 * u.rad,
            u0 * u.m / u.s,
            u1 * u.rad / u.s,
            u2 * u.rad / u.s,
        )

    def bl_differential(self, **kwargs):
        """
        Converts to Boyer-Lindquist Coordinates

        Parameters
        ----------
        **kwargs : dict
            Keyword Arguments
            Expects two arguments, ``M and ``a``, as described below

        Other Parameters
        ----------------
        M : float
            Mass of the gravitating body, \
            around which, spacetime has been defined
        a : float
            Spin Parameter of the gravitating body, \
            around which, spacetime has been defined

        Returns
        -------
        ~einsteinpy.coordinates.differential.BoyerLindquistDifferential
            Boyer-Lindquist representation of velocity

        """
        M, a = kwargs["M"], kwargs["a"]
        e0, e1, e2, e3, u0, u1, u2 = self.convert_bl(M=M, a=a)
        return BoyerLindquistDifferential(
            e0 * u.s,
            e1 * u.m,
            e2 * u.rad,
            e3 * u.rad,
            u0 * u.m / u.s,
            u1 * u.rad / u.s,
            u2 * u.rad / u.s,
        )


class SphericalDifferential(BaseDifferential, SphericalConversion):
    """
    Class for defining 3-Velocity & 4-Velocity in Spherical Polar Coordinates \
    using SI units

    """

    @u.quantity_input(
        e0=u.s,
        e1=u.m,
        e2=u.rad,
        e3=u.rad,
        u0=u.m / u.s,
        u1=u.rad / u.s,
        u2=u.rad / u.s,
    )
    def __init__(self, e0, e1, e2, e3, u0, u1, u2):
        super().__init__(e0, e1, e2, e3, u0, u1, u2)
        self.system = "Spherical"

    @property
    def v_t(self):
        """
        Returns the Timelike component of 4-Velocity

        """
        return self._v_t

    @v_t.setter
    def v_t(self, args):
        """
        Sets the value of the Time-like component of 4-Velocity

        Parameters
        ----------
        args : tuple
            1-tuple containing the ~einsteinpy.metric.* object, \
            in which the coordinates are defined

        Raises
        ------
        CoordinateError
            If ``metric`` object has been instantiated with a coordinate system, \
            other than Sperical Polar Coordinates.

        """
        g = args[0]
        if self.system != g.coords.system:
            raise CoordinateError(
                f"Metric object has been instantiated with a coordinate system, ( {g.coords.system} )"
                " other than Spherical Polar Coordinates."
            )

        g_cov_mat = g.metric_covariant(self.position())

        v_t = v0(g_cov_mat, self.u0.si.value, self.u1.si.value, self.u2.si.value)

        self._v_t = v_t * u.m / u.s

    def cartesian_differential(self, **kwargs):
        """
        Converts to Cartesian Coordinates

        Parameters
        ----------
        **kwargs : dict
            Keyword Arguments

        Returns
        -------
        ~einsteinpy.coordinates.differential.CartesianDifferential
            Cartesian representation of velocity

        """
        e0, e1, e2, e3, u0, u1, u2 = self.convert_cartesian()
        return CartesianDifferential(
            e0 * u.s,
            e1 * u.m,
            e2 * u.m,
            e3 * u.m,
            u0 * u.m / u.s,
            u1 * u.m / u.s,
            u2 * u.m / u.s,
        )

    def bl_differential(self, **kwargs):
        """
        Converts to Boyer-Lindquist coordinates

        Parameters
        ----------
        **kwargs : dict
            Keyword Arguments
            Expects two arguments, ``M and ``a``, as described below

        Other Parameters
        ----------------
        M : float
            Mass of the gravitating body, \
            around which, spacetime has been defined
        a : float
            Spin Parameter of the gravitating body, \
            around which, spacetime has been defined

        Returns
        -------
        ~einsteinpy.coordinates.differential.BoyerLindquistDifferential
            Boyer-Lindquist representation of velocity

        """
        M, a = kwargs["M"], kwargs["a"]
        e0, e1, e2, e3, u0, u1, u2 = self.convert_bl(M=M, a=a)
        return BoyerLindquistDifferential(
            e0 * u.s,
            e1 * u.m,
            e2 * u.rad,
            e3 * u.rad,
            u0 * u.m / u.s,
            u1 * u.rad / u.s,
            u2 * u.rad / u.s,
        )


class BoyerLindquistDifferential(BaseDifferential, BoyerLindquistConversion):
    """
    Class for defining 3-Velocity & 4-Velocity in Boyer-Lindquist Coordinates \
    using SI units

    """

    @u.quantity_input(
        e0=u.s,
        e1=u.m,
        e2=u.rad,
        e3=u.rad,
        u0=u.m / u.s,
        u1=u.rad / u.s,
        u2=u.rad / u.s,
    )
    def __init__(self, e0, e1, e2, e3, u0, u1, u2):
        super().__init__(e0, e1, e2, e3, u0, u1, u2)
        self.system = "BoyerLindquist"

    @property
    def v_t(self):
        """
        Returns the Timelike component of 4-Velocity

        """
        return self._v_t

    @v_t.setter
    def v_t(self, args):
        """
        Sets the value of the Time-like component of 4-Velocity

        Parameters
        ----------
        args : tuple
            1-tuple containing the ~einsteinpy.metric.* object, \
            in which the coordinates are defined

        Raises
        ------
        CoordinateError
            If ``metric`` object has been instantiated with a coordinate system, \
            other than Boyer-Lindquist Coordinates.

        """
        g = args[0]
        if self.system != g.coords.system:
            raise CoordinateError(
                "Metric object has been instantiated with a coordinate system, ( {g.coords.system} )"
                " other than Boyer-Lindquist Coordinates."
            )

        g_cov_mat = g.metric_covariant(self.position())

        v_t = v0(g_cov_mat, self.u0.si.value, self.u1.si.value, self.u2.si.value)

        self._v_t = v_t * u.m / u.s

    def cartesian_differential(self, **kwargs):
        """
        Converts to Cartesian Coordinates

        Parameters
        ----------
        **kwargs : dict
            Keyword Arguments
            Expects two arguments, ``M and ``a``, as described below

        Other Parameters
        ----------------
        M : float
            Mass of the gravitating body, \
            around which, spacetime has been defined
        a : float
            Spin Parameter of the gravitating body, \
            around which, spacetime has been defined

        Returns
        -------
        ~einsteinpy.coordinates.differentia.CartesianDifferential
            Cartesian representation of velocity

        """
        M, a = kwargs["M"], kwargs["a"]
        e0, e1, e2, e3, u0, u1, u2 = self.convert_cartesian(M=M, a=a)
        return CartesianDifferential(
            e0 * u.s,
            e1 * u.m,
            e2 * u.m,
            e3 * u.m,
            u0 * u.m / u.s,
            u1 * u.m / u.s,
            u2 * u.m / u.s,
        )

    def spherical_differential(self, **kwargs):
        """
        Converts to Spherical Polar Coordinates

        Parameters
        ----------
        **kwargs : dict
            Keyword Arguments
            Expects two arguments, ``M and ``a``, as described below

        Other Parameters
        ----------------
        M : float
            Mass of the gravitating body, \
            around which, spacetime has been defined
        a : float
            Spin Parameter of the gravitating body, \
            around which, spacetime has been defined

        Returns
        -------
        ~einsteinpy.coordinates.differentia.SphericalDifferential
            Spherical representation of velocity

        """
        M, a = kwargs["M"], kwargs["a"]
        e0, e1, e2, e3, u0, u1, u2 = self.convert_spherical(M=M, a=a)
        return SphericalDifferential(
            e0 * u.s,
            e1 * u.m,
            e2 * u.rad,
            e3 * u.rad,
            u0 * u.m / u.s,
            u1 * u.rad / u.s,
            u2 * u.rad / u.s,
        )
