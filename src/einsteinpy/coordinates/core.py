import numpy as np

from einsteinpy.coordinates.conversion import (
    BoyerLindquistConversion,
    CartesianConversion,
    SphericalConversion,
)
from einsteinpy.coordinates.utils import four_position


class Cartesian(CartesianConversion):
    """
    Class for defining 3-Position & 4-Position in Cartesian Coordinates \
    using SI units

    """

    def __init__(self, t, x, y, z):
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

        """
        self.t = t
        self.x = x
        self.y = y
        self.z = z
        super().__init__(t, x, y, z)
        self.system = "Cartesian"
        self._dimension = {
            "t": self.t,
            "x": self.x,
            "y": self.y,
            "z": self.z,
            "system": self.system,
        }
        self._dimension_order = ("t", "x", "y", "z")

    def __str__(self):
        return f"Cartesian Coordinates: \n \
            t = ({self.t}), x = ({self.x}), y = ({self.y}), z = ({self.z})"

    def __repr__(self):
        return f"Cartesian Coordinates: \n \
            t = ({self.t}), x = ({self.x}), y = ({self.y}), z = ({self.z})"

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
        ~numpy.ndarray :
            Array, containing Position 4-Vector in SI units

        """
        x4 = four_position(self.t, self.x, self.y, self.z)

        return x4

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
        t, r, theta, phi = self.convert_spherical()

        return Spherical(t, r, theta, phi)

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
        t, r, theta, phi = self.convert_bl(M=M, a=a)

        return BoyerLindquist(t, r, theta, phi)


class Spherical(SphericalConversion):
    """
    Class for defining 3-Position & 4-Position in Spherical Polar Coordinates \
    using SI units

    """

    def __init__(self, t, r, theta, phi):
        """
        Constructor

        Parameters
        ----------
        t : float
            Time
        r : float
            r-Component of 3-Position
        theta : float
            theta-Component of 3-Position
        phi : float
            phi-Component of 3-Position

        """
        self.t = t
        self.r = r
        self.theta = theta
        self.phi = phi
        super().__init__(t, r, theta, phi)
        self.system = "Spherical"
        self._dimension = {
            "t": self.t,
            "r": self.r,
            "theta": self.theta,
            "phi": self.phi,
            "system": self.system,
        }
        self._dimension_order = ("t", "r", "theta", "phi")

    def __str__(self):
        return f"Spherical Polar Coordinates: \n \
            t = ({self.t}), r = ({self.r}), theta = ({self.theta}), phi = ({self.phi})"

    def __repr__(self):
        return f"Spherical Polar Coordinates: \n \
            t = ({self.t}), r = ({self.r}), theta = ({self.theta}), phi = ({self.phi})"

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
        ~numpy.ndarray :
            Array, containing Position 4-Vector in SI units

        """
        x4 = four_position(self.t, self.r, self.theta, self.phi)

        return x4

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
        t, x, y, z = self.convert_cartesian()
        return Cartesian(t, x, y, z)

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
        t, r, theta, phi = self.convert_bl(M=M, a=a)
        return BoyerLindquist(t, r, theta, phi)


class BoyerLindquist(BoyerLindquistConversion):
    """
    Class for defining 3-Position & 4-Position in Boyer-Lindquist Coordinates \
    using SI units

    """

    def __init__(self, t, r, theta, phi):
        """
        Constructor

        Parameters
        ----------
        t : float
            Time
        r : float
            r-Component of 3-Position
        theta : float
            theta-Component of 3-Position
        phi : float
            phi-Component of 3-Position

        """
        self.t = t
        self.r = r
        self.theta = theta
        self.phi = phi
        super().__init__(t, r, theta, phi)
        self.system = "BoyerLindquist"
        self._dimension = {
            "t": self.t,
            "r": self.r,
            "theta": self.theta,
            "phi": self.phi,
            "system": self.system,
        }
        self._dimension_order = ("t", "r", "theta", "phi")

    def __str__(self):
        return f"Boyer-Lindquist Coordinates: \n \
            t = ({self.t}), r = ({self.r}), theta = ({self.theta}), phi = ({self.phi})"

    def __repr__(self):
        return f"Boyer-Lindquist Coordinates: \n \
            t = ({self.t}), r = ({self.r}), theta = ({self.theta}), phi = ({self.phi})"

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
        ~numpy.ndarray :
            Array, containing Position 4-Vector in SI units

        """
        x4 = four_position(self.t, self.r, self.theta, self.phi)

        return x4

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
        t, x, y, z = self.convert_cartesian(M=M, a=a)
        return Cartesian(t, x, y, z)

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
        t, r, theta, phi = self.convert_spherical(M=M, a=a)
        return Spherical(t, r, theta, phi)
