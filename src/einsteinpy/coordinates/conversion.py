from einsteinpy.coordinates.utils import (
    bl_to_cartesian_fast,
    cartesian_to_bl_fast,
    cartesian_to_spherical_fast,
    spherical_to_cartesian_fast,
)
from einsteinpy.metric import BaseMetric
from typing import Union, Tuple, Optional, Dict


class CoordinateConversionBase:
    """
    Base class for conversion to and from different coordinate systems in SI units
    """

    def __init__(
        self,
        t: float,
        x: float,
        y: float,
        z: float,
        v_x: Optional[float] = None,
        v_y: Optional[float] = None,
        v_z: Optional[float] = None,
    ):
        """
        Constructor

        Parameters
        ----------
        t : float
            Time
        x : float
            x-Component of 3-Position
        y : float
            y-Component of 3-Position
        z : float
            z-Component of 3-Position
        v_x : float, optional
            x-Component of 3-Velocity
        v_y : float, optional
            y-Component of 3-Velocity
        v_z : float, optional
            z-Component of 3-Velocity
        """
        self.t_si = t
        self.x_si = x
        self.y_si = y
        self.z_si = z
        self.v_x_si = v_x
        self.v_y_si = v_y
        self.v_z_si = v_z
        self._velocities_provided = not (
            (v_x is None) or (v_y is None) or (v_z is None)
        )

    def values(self) -> Union[Tuple[float, float, float, float], Tuple[float, float, float, float, float, float, float]]:
        """
        Returns components of the coordinates in SI units

        Returns
        -------
        tuple
            4-Tuple, containing ``t, x, y, z`` in SI units
            or 7-tuple, containing ``t, x, y, z, v_x, v_y, v_z`` \
            in SI units
        """
        if self._velocities_provided:
            return (
                self.t_si,
                self.x_si,
                self.y_si,
                self.z_si,
                self.v_x_si,
                self.v_y_si,
                self.v_z_si,
            )

        return self.t_si, self.x_si, self.y_si, self.z_si


class CartesianConversion(CoordinateConversionBase):
    """
    Class for conversion to and from Cartesian Coordinates in SI units
    """

    def convert_spherical(self) -> Tuple[float, float, float, float]:
        """
        Converts to Spherical Polar Coordinates

        Returns
        -------
        tuple
            4-Tuple containing the components in Spherical Polar Coordinates
        """
        return cartesian_to_spherical_fast(*self.values())

    def convert_bl(self, M: float, a: float) -> Tuple[float, float, float, float]:
        """
        Converts to Boyer-Lindquist Coordinates

        Parameters
        ----------
        M : float
            Mass of the gravitating body
        a : float
            Spin Parameter of the gravitating body

        Returns
        -------
        tuple
            4-Tuple containing the components in Boyer-Lindquist Coordinates
        """
        alpha = BaseMetric.alpha(M=M, a=a)
        return cartesian_to_bl_fast(*self.values(), alpha)


class SphericalConversion(CoordinateConversionBase):
    """
    Class for conversion to and from Spherical Polar Coordinates in SI units
    """

    def convert_cartesian(self) -> Tuple[float, float, float, float]:
        """
        Converts to Cartesian Coordinates

        Returns
        -------
        tuple
            4-Tuple containing the components in Cartesian Coordinates
        """
        return spherical_to_cartesian_fast(*self.values())

    def convert_bl(self, M: float, a: float) -> Tuple[float, float, float, float]:
        """
        Converts to Boyer-Lindquist Coordinates

        Parameters
        ----------
        M : float
            Mass of the gravitating body
        a : float
            Spin Parameter of the gravitating body

        Returns
        -------
        tuple
            4-Tuple containing the components in Boyer-Lindquist Coordinates
        """
        transformed_cartesian = self.convert_cartesian()
        cart = CartesianConversion(*transformed_cartesian)
        return cart.convert_bl(M=M, a=a)


class BoyerLindquistConversion(CoordinateConversionBase):
    """
    Class for conversion to and from Boyer-Lindquist Coordinates in SI units
    """

    def convert_cartesian(self, M: float, a: float) -> Tuple[float, float, float, float]:
        """
        Converts to Cartesian Coordinates

        Parameters
        ----------
        M : float
            Mass of the gravitating body
        a : float
            Spin Parameter of the gravitating body

        Returns
        -------
        tuple
            4-Tuple containing the components in Cartesian Coordinates
        """
        alpha = BaseMetric.alpha(M=M, a=a)
        return bl_to_cartesian_fast(*self.values(), alpha)

    def convert_spherical(self, M: float, a: float) -> Tuple[float, float, float, float]:
        """
        Converts to Spherical Polar Coordinates

        Parameters
        ----------
        M : float
            Mass of the gravitating body
        a : float
            Spin Parameter of the gravitating body

        Returns
        -------
        tuple
            4-Tuple containing the components in Spherical Polar Coordinates
        """
        transformed_cartesian = self.convert_cartesian(M=M, a=a)
        cart = CartesianConversion(*transformed_cartesian)
        return cart.convert_spherical()
