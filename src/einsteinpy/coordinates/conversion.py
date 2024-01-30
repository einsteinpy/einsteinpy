from einsteinpy.coordinates.utils import (
    bl_to_cartesian_fast,
    cartesian_to_bl_fast,
    cartesian_to_spherical_fast,
    spherical_to_cartesian_fast,
)
from einsteinpy.metric import BaseMetric


class CartesianConversion:
    """
    Class for conversion to and from Cartesian Coordinates in SI units

    """

    def __init__(self, e0, e1, e2, e3, u0=None, u1=None, u2=None):
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
        u0 : float, optional
            x-Component of 3-Velocity
        u1 : float, optional
            y-Component of 3-Velocity
        u2 : float, optional
            z-Component of 3-Velocity

        """
        self.e0_si = e0
        self.e1_si = e1
        self.e2_si = e2
        self.e3_si = e3
        self.u0_si = u0
        self.u1_si = u1
        self.u2_si = u2
        self._velocities_provided = not ((u0 is None) or (u1 is None) or (u2 is None))

    def values(self):
        """
        Returns components of the coordinates in SI units

        Returns
        -------
        tuple
            4-Tuple, containing ``e0, e1, e2, e3`` in SI units
            or 7-tuple, containing ``e0, e1, e2, e3, u0, u1, u2`` \
            in SI units

        """
        if self._velocities_provided:
            return (
                self.e0_si,
                self.e1_si,
                self.e2_si,
                self.e3_si,
                self.u0_si,
                self.u1_si,
                self.u2_si,
            )

        return self.e0_si, self.e1_si, self.e2_si, self.e3_si

    def convert_spherical(self, **kwargs):
        """
        Converts to Spherical Polar Coordinates

        Other Parameters
        ----------------
        **kwargs : dict
            Keyword Arguments

        Returns
        -------
        tuple
            4-Tuple or 7-Tuple, containing the components in \
            Spherical Polar Coordinates

        """
        return cartesian_to_spherical_fast(
            self.e0_si,
            self.e1_si,
            self.e2_si,
            self.e3_si,
            self.u0_si,
            self.u1_si,
            self.u2_si,
            self._velocities_provided,
        )

    def convert_bl(self, **kwargs):
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
        tuple
            4-Tuple or 7-Tuple, containing the components in \
            Boyer-Lindquist Coordinates

        Raises
        ------
        KeyError
            If ``kwargs`` does not contain both ``M`` \
            and ``a`` as keyword arguments

        """
        try:
            M, a = kwargs["M"], kwargs["a"]
        except KeyError:
            raise KeyError(
                "Two keyword arguments are expected: Mass, 'M' and Spin Parameter, 'a'."
            )

        alpha = BaseMetric.alpha(M=M, a=a)

        return cartesian_to_bl_fast(
            self.e0_si,
            self.e1_si,
            self.e2_si,
            self.e3_si,
            alpha,
            self.u0_si,
            self.u1_si,
            self.u2_si,
            self._velocities_provided,
        )


class SphericalConversion:
    """
    Class for conversion to and from Spherical Polar Coordinates in SI units

    """

    def __init__(self, e0, e1, e2, e3, u0=None, u1=None, u2=None):
        """
        Constructor

        Parameters
        ----------
        e0 : float
            Time
        e1 : float
            r-Component of 3-Position
        e2 : float
            theta-Component of 3-Position
        e3 : float
            phi-Component of 3-Position
        u0 : float, optional
            r-Component of 3-Velocity
        u1 : float, optional
            theta-Component of 3-Velocity
        u2 : float, optional
            phi-Component of 3-Velocity

        """
        self.e0_si = e0
        self.e1_si = e1
        self.e2_si = e2
        self.e3_si = e3
        self.u0_si = u0
        self.u1_si = u1
        self.u2_si = u2
        self._velocities_provided = not ((u0 is None) or (u1 is None) or (u2 is None))

    def values(self):
        """
        Returns components of the coordinates

        Returns
        -------
        tuple
            4-Tuple containing ``t, r, theta, phi`` in SI units
            or 7-tuple, containing ``t, r, theta, phi, v_r, v_th, v_p`` \
            in SI units

        """
        if self._velocities_provided:
            return (
                self.e0_si,
                self.e1_si,
                self.e2_si,
                self.e3_si,
                self.u0_si,
                self.u1_si,
                self.u2_si,
            )

        return self.e0_si, self.e1_si, self.e2_si, self.e3_si

    def convert_cartesian(self, **kwargs):
        """
        Converts to Cartesian Coordinates

        Other Parameters
        ----------------
        **kwargs : dict
            Keyword Arguments

        Returns
        -------
        tuple
            4-Tuple or 7-Tuple, containing the components in \
            Cartesian Coordinates

        """
        return spherical_to_cartesian_fast(
            self.e0_si,
            self.e1_si,
            self.e2_si,
            self.e3_si,
            self.u0_si,
            self.u1_si,
            self.u2_si,
            self._velocities_provided,
        )

    def convert_bl(self, **kwargs):
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
        tuple
            4-Tuple or 7-Tuple, containing the components in \
            Boyer-Lindquist Coordinates

        Raises
        ------
        KeyError
            If ``kwargs`` does not contain both ``M`` \
            and ``a`` as keyword arguments

        """
        try:
            M, a = kwargs["M"], kwargs["a"]
        except KeyError:
            raise KeyError(
                "Two keyword arguments are expected: Mass, 'M' and Spin Parameter, 'a'."
            )

        transformed_cartesian = self.convert_cartesian()
        cart = CartesianConversion(*transformed_cartesian)

        return cart.convert_bl(M=M, a=a)


class BoyerLindquistConversion:
    """
    Class for conversion to and from Boyer-Lindquist Coordinates in SI units

    """

    def __init__(self, e0, e1, e2, e3, u0=None, u1=None, u2=None):
        """
        Constructor

        Parameters
        ----------
        e0 : float
            Time
        e1 : float
            r-Component of 3-Position
        e2 : float
            theta-Component of 3-Position
        e3 : float
            phi-Component of 3-Position
        u0 : float, optional
            r-Component of 3-Velocity
        u1 : float, optional
            theta-Component of 3-Velocity
        u2 : float, optional
            phi-Component of 3-Velocity

        """
        self.e0_si = e0
        self.e1_si = e1
        self.e2_si = e2
        self.e3_si = e3
        self.u0_si = u0
        self.u1_si = u1
        self.u2_si = u2
        self._velocities_provided = not ((u0 is None) or (u1 is None) or (u2 is None))

    def values(self):
        """
        Returns components of the coordinates

        Returns
        -------
        tuple
            4-Tuple containing ``t, r, theta, phi`` in SI units
            or 7-tuple, containing ``t, r, theta, phi, v_r, v_th, v_p`` \
            in SI units

        """
        if self._velocities_provided:
            return (
                self.e0_si,
                self.e1_si,
                self.e2_si,
                self.e3_si,
                self.u0_si,
                self.u1_si,
                self.u2_si,
            )

        return self.t_si, self.r_si, self.th_si, self.p_si

    def convert_cartesian(self, **kwargs):
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
        tuple
            4-Tuple or 7-Tuple, containing the components in \
            Cartesian Coordinates

        Raises
        ------
        KeyError
            If ``kwargs`` does not contain both ``M`` \
            and ``a`` as keyword arguments

        """
        try:
            M, a = kwargs["M"], kwargs["a"]
        except KeyError:
            raise KeyError(
                "Two keyword arguments are expected: Mass, 'M' and Spin Parameter, 'a'."
            )

        alpha = BaseMetric.alpha(M=M, a=a)

        return bl_to_cartesian_fast(
            self.e0_si,
            self.e1_si,
            self.e2_si,
            self.e3_si,
            alpha,
            self.u0_si,
            self.u1_si,
            self.u2_si,
            self._velocities_provided,
        )

    def convert_spherical(self, **kwargs):
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
        tuple
            4-Tuple or 7-Tuple, containing the components in \
            Spherical Polar Coordinates

        Raises
        ------
        KeyError
            If ``kwargs`` does not contain both ``M`` \
            and ``a`` as keyword arguments

        """
        try:
            M, a = kwargs["M"], kwargs["a"]
        except KeyError:
            raise KeyError(
                "Two keyword arguments are expected: Mass, 'M' and Spin Parameter, 'a'."
            )

        transformed_cartesian = self.convert_cartesian(M=M, a=a)
        cart = CartesianConversion(*transformed_cartesian)

        return cart.convert_spherical()
