import numpy as np
from astropy import units as u

from einsteinpy import constant
from einsteinpy.coordinates.conversion import (
    BoyerLindquistConversion,
    CartesianConversion,
    SphericalConversion,
)

_c = constant.c.value


class Cartesian(CartesianConversion):
    """
    Class for defining 3-Position & 4-Position in Cartesian Coordinates \
    using SI units

    """

    @u.quantity_input(t=u.s, x=u.m, y=u.m, z=u.m)
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
        super().__init__(t.si.value, x.si.value, y.si.value, z.si.value)
        self.t = t
        self.x = x
        self.y = y
        self.z = z
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
        tuple
            4-Tuple, containing Position 4-Vector in SI units

        """
        return (_c * self.t.si.value, self.x.si.value, self.y.si.value, self.z.si.value)

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

        return Spherical(t * u.s, r * u.m, theta * u.rad, phi * u.rad)

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

        return BoyerLindquist(t * u.s, r * u.m, theta * u.rad, phi * u.rad)


class Spherical(SphericalConversion):
    """
    Class for defining 3-Position & 4-Position in Spherical Polar Coordinates \
    using SI units

    """

    @u.quantity_input(t=u.s, r=u.m, theta=u.rad, phi=u.rad)
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
        super().__init__(t.si.value, r.si.value, theta.si.value, phi.si.value)
        self.t = t
        self.r = r
        self.theta = theta
        self.phi = phi
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
        tuple :
            4-Tuple, containing Position 4-Vector in SI units

        """
        return (
            _c * self.t.si.value,
            self.r.si.value,
            self.theta.si.value,
            self.phi.si.value,
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
        t, x, y, z = self.convert_cartesian()

        return Cartesian(t * u.s, x * u.m, y * u.m, z * u.m)

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

        return BoyerLindquist(t * u.s, r * u.m, theta * u.rad, phi * u.rad)

    def to_eddingtonfinkelstein(self, **kwargs):
        """
        Method for conversion to ingoing Eddington-Finkelstein Coordinates

        Parameters
        ---------
        **kwargs : dict
            Keyword Arguments
            Expects one argument, ``M``, as described below

        Other Parameters
        ---------------
        M : float
            Mass of gravitating body
            Required to calculate r_star, the Tortoise coordinate \
            to transform t to v

        Returns
        ------
        ~einsteinpy.coordinates.core.IngoingEddingtonFinkelstein
            Ingoing Eddington-Finkelstein representation of the Spherical Polar Coordinates

        """
        M = kwargs["M"]
        v, r, theta, phi = self.convert_ief(M=M)

        return IngoingEddingtonFinkelstein(v * u.s, r * u.rad, theta * u.rad, phi * u.rad)


class BoyerLindquist(BoyerLindquistConversion):
    """
    Class for defining 3-Position & 4-Position in Boyer-Lindquist Coordinates \
    using SI units

    """

    @u.quantity_input(t=u.s, r=u.m, theta=u.rad, phi=u.rad)
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
        super().__init__(t.si.value, r.si.value, theta.si.value, phi.si.value)
        self.t = t
        self.r = r
        self.theta = theta
        self.phi = phi
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
        tuple :
            4-Tuple, containing Position 4-Vector in SI units

        """
        return (
            _c * self.t.si.value,
            self.r.si.value,
            self.theta.si.value,
            self.phi.si.value,
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
        t, x, y, z = self.convert_cartesian(M=M, a=a)

        return Cartesian(t * u.s, x * u.m, y * u.m, z * u.m)

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

        return Spherical(t * u.s, r * u.m, theta * u.rad, phi * u.rad)


class KerrSchild(KerrSchildConversion):
    """
    Class for defining 3-position and 4-position in Kerr-Schild coordinates \
    using SI units

    """
    @u.quantity_input(t=u.s, x=u.m, y=u.m, z=u.m)
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
        super().__init__(t.si.value, x.si.value, y.si.value, z.si.value)
        self.t = t
        self.x = x
        self.y = y
        self.z = z
        self.system = "Kerr-Schild"
        self._dimension = {
            "t": self.t,
            "x": self.x,
            "y": self.y,
            "z": self.z,
            "system": self.system,
        }
        self._dimension_order = ("t", "x", "y", "z")

    def __str__(self):
        return f"Kerr-Schild Coordinates: \n \
            t = ({self.t}), x = ({self.x}), y = ({self.y}), z = ({self.z})"

    def __repr__(self):
        return f"Kerr-Schild Coordinates: \n \
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
        tuple :
            4-Tuple, containing Position 4-Vector in SI units

        """
        return (
            _c * self.t.si.value,
            self.x.si.value,
            self.y.si.value,
            self.z.si.value,
        )

    def to_eddingtonfinkelstein(self, **kwargs):
        """
        Method for conversion to ingoing Eddington-Finkelstein Coordinates

        Parameters
        ---------
        **kwargs : dict
            Keyword Arguments
            Expects two arguments, ``M`` and ``a``, as described below

        Other Parameters
        ---------------
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
        ------
        ~einsteinpy.coordinates.core.IngoingEddingtonFinkelstein
            Ingoing Eddington-Finkelstein representation of the Kerr-Schild Coordinates

        """
        M, a = kwargs["M"], kwargs["a"]
        v, r, theta, phi = self.convert_ief(M=M, a=a)

        return IngoingEddingtonFinkelstein(v * u.s, r * u.rad, theta * u.rad, phi * u.rad)


class IngoingEddingtonFinkelstein(IngoingEddingtonFinkelsteinConversion):
    """
    Class for defining 3-position and 4-position in Eddington-Finkelstein coordinates \
    using SI units

    """
    @u.quantity_input(v=u.s, r=u.m, theta=u.rad, phi=u.rad)
    def __init__(self, v, r, theta, phi):
        """
        Constructor

        Parameters
        ----------
        v : float
            t + (r_star/c) where t is Time, r_star is Tortoise Coordinate
        r : float
            r-Component of 3-Position
        theta : float
            theta-Component of 3-Position
        phi : float
            phi-Component of 3-Position

        """
        super().__init__(v.si.value, r.si.value, theta.si.value, phi.si.value)
        self.v = v
        self.r = r
        self.theta = theta
        self.phi = phi
        self.system = "Ingoing Eddington-Finkelstein"
        self._dimension = {
            "v": self.v,
            "r": self.r,
            "theta": self.theta,
            "phi": self.phi,
            "system": self.system,
        }
        self._dimension_order = ("v", "r", "theta", "phi")

    def __str__(self):
        return f"Ingoing Eddington-Finkelstein Coordinates: \n \
            v = ({self.v}), r = ({self.r}), theta = ({self.theta}), phi = ({self.phi})"

    def __repr__(self):
        return f"Ingoing Eddington-Finkelstein Coordinates: \n \
            v = ({self.v}), r = ({self.r}), theta = ({self.theta}), phi = ({self.phi})"

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
            _c * self.v.si.value,
            self.r.si.value,
            self.theta.si.value,
            self.phi.si.value,
        )

    def to_spherical(self, **kwargs):
        """
        Method for conversion to Spherical Polar Coordinates

        Parameters
        ----------
        **kwargs : dict
            Keyword Arguments
            Expects one argument, ``M``, as described below

        Other Parameters
        ----------------
        M : float
            Mass of gravitating body
            Required to calculate r_star, the Tortoise coordinate \
            to transform v to t

        Returns
        -------
        ~einsteinpy.coordinates.core.Spherical
            Spherical representation of the ingoing Eddington-Finkelstein \
            Coordinates

        """
        M = kwargs["M"]
        t, r, theta, phi = self.convert_spherical(M=M)

        return Spherical(t * u.s, r * u.m, theta * u.rad, phi * u.rad)

    def to_kerrschild(self, **kwargs):
        """
        Method for conversion to Kerr-Schild Coordinates

        Parameters
        ----------
        **kwargs : dict
            Keyword Arguments
            Expects two argument, ``M`` and ``a``, as described below

        Other Parameters
        ---------------
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
        ------
        ~einsteinpy.coordinates.core.KerrSchild
            Kerr-Schild representation of the ingoing Eddington-Finkelstein Coordinates

        """
        M, a = kwargs["M"], kwargs["a"]
        t, x, y, z = self.convert_ks(M=M, a=a)

        return KerrSchild(t * u.s, x * u.m, y * u.m, z * u.m)