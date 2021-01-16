import astropy.units as u
import numpy as np

from einsteinpy.coordinates.conversion import (
    BoyerLindquistConversion,
    CartesianConversion,
    SphericalConversion,
)


class Cartesian(CartesianConversion):
    """
    Class for Cartesian Coordinates and related transformations.
    """

    @u.quantity_input(x=u.km, y=u.km, z=u.km)
    def __init__(self, x, y, z):
        """
        Constructor.

        Parameters
        ----------
        x : ~astropy.units.quantity.Quantity
        y : ~astropy.units.quantity.Quantity
        z : ~astropy.units.quantity.Qauntity

        """
        self.x = x
        self.y = y
        self.z = z
        super().__init__(x.si.value, y.si.value, z.si.value)
        self.system = "Cartesian"
        self._dimension = {"x": self.x, "y": self.y, "z": self.z, "system": self.system}
        self._dimension_order = ("x", "y", "z")

    def __repr__(self):
        return "Cartesian x: {}, y: {}, z: {}".format(self.x, self.y, self.z)

    def __str__(self):
        return self.__repr__()

    def __getitem__(self, item):
        """
        Method to return coordinates.
        Objects are subsctiptable with both explicit names of parameters
        and integer indices.

        Parameters
        ----------
        item : str or int
                Name of the parameter or its index.
                If ``'system'``, Name of coordinate is returned.
        """
        if isinstance(item, (int, np.integer)):
            return self._dimension[self._dimension_order[item]]
        return self._dimension[item]

    def si_values(self):
        """
        Function for returning values in SI units.

        Returns
        -------
        ~numpy.ndarray
            Array containing values in SI units (m, m, m)

        """
        element_list = [self.x.to(u.m), self.y.to(u.m), self.z.to(u.m)]
        return np.array([e.value for e in element_list], dtype=float)

    def norm(self):
        """
        Function for finding euclidean norm of a vector.

        Returns
        -------
        ~astropy.units.quantity.Quantity
            Euclidean norm with units.

        """
        return np.sqrt(self.x ** 2 + self.y ** 2 + self.z ** 2)

    def dot(self, target):
        """
        Dot product of two vectors.

        Parameters
        ----------
        target: ~einsteipy.coordinates.core.Cartesian

        Returns
        -------
        ~astropy.units.quantity.Quantity
            Dot product with units

        """

        x = self.x * target.x
        y = self.y * target.y
        z = self.z * target.z

        return x + y + z

    def to_spherical(self):
        """
        Method for conversion to spherical coordinates.

        Returns
        -------
        ~einsteinpy.coordinates.core.Spherical
            Spherical representation of the Cartesian Coordinates.

        """
        r, theta, phi = self.convert_spherical()
        return Spherical(r * u.m, theta * u.rad, phi * u.rad)

    @u.quantity_input(a=u.km)
    def to_bl(self, a):
        """
        Method for conversion to boyer-lindquist coordinates.

        Parameters
        ----------
        a : ~astropy.units.quantity.Quantity
            a = J/Mc , the angular momentum per unit mass of the black hole per speed of light.

        Returns
        -------
        ~einsteinpy.coordinates.core.BoyerLindquist
            BL representation of the Cartesian Coordinates.

        """
        r, theta, phi, a = self.convert_bl(a.si.value)
        return BoyerLindquist(r * u.m, theta * u.rad, phi * u.rad, a * u.m)


class Spherical(SphericalConversion):
    """
    Class for Spherical Coordinates and related transformations.
    """

    @u.quantity_input(r=u.km, theta=u.rad, phi=u.rad)
    def __init__(self, r, theta, phi):
        """
        Constructor.

        Parameters
        ----------
        r : ~astropy.units.quantity.Quantity
        theta : ~astropy.units.quantity.Quantity
        phi : ~astropy.units.quantity.Quantity

        """
        self.r = r
        self.theta = theta
        self.phi = phi
        super().__init__(r.si.value, theta.si.value, phi.si.value)
        self.system = "Spherical"
        self._dimension = {
            "r": self.r,
            "theta": self.theta,
            "phi": self.phi,
            "system": self.system,
        }
        self._dimension_order = ("r", "theta", "phi")

    def __repr__(self):
        return "Spherical r: {}, theta: {}, phi: {}".format(
            self.r, self.theta, self.phi
        )

    def __str__(self):
        return self.__repr__()

    def __getitem__(self, item):
        """
        Method to return coordinates.
        Objects are subsctiptable with both explicit names of parameters
        and integer indices.

        Parameters
        ----------
        item : str or int
                Name of the parameter or its index.
                If ``'system'``, Name of coordinate is returned.
        """
        if isinstance(item, (int, np.integer)):
            return self._dimension[self._dimension_order[item]]
        return self._dimension[item]

    def si_values(self):
        """
        Function for returning values in SI units.

        Returns
        -------
        ~numpy.ndarray
            Array containing values in SI units (m, rad, rad)

        """
        element_list = [self.r.to(u.m), self.theta.to(u.rad), self.phi.to(u.rad)]
        return np.array([e.value for e in element_list], dtype=float)

    def to_cartesian(self):
        """
        Method for conversion to cartesian coordinates.

        Returns
        -------
        ~einsteinpy.coordinates.core.Cartesian
            Cartesian representation of the Spherical Coordinates.

        """
        x, y, z = self.convert_cartesian()
        return Cartesian(x * u.m, y * u.m, z * u.m)

    @u.quantity_input(a=u.km)
    def to_bl(self, a):
        """
        Method for conversion to boyer-lindquist coordinates.

        Parameters
        ----------
        a : ~astropy.units.quantity.Quantity
            a = J/Mc , the angular momentum per unit mass of the black hole per speed of light.

        Returns
        -------
        ~einsteinpy.coordinates.core.BoyerLindquist
            BL representation of the Spherical Coordinates.

        """
        r, theta, phi, a = self.convert_bl(a.si.value)
        return BoyerLindquist(r * u.m, theta * u.rad, phi * u.rad, a * u.m)


class BoyerLindquist(BoyerLindquistConversion):
    """
    Class for Spherical Coordinates and related transformations.
    """

    @u.quantity_input(r=u.km, theta=u.rad, phi=u.rad, a=u.km)
    def __init__(self, r, theta, phi, a):
        """
        Constructor.

        Parameters
        ----------
        r : ~astropy.units.quantity.Quantity
        theta : ~astropy.units.quantity.Quantity
        phi : ~astropy.units.quantity.Quantity
        a : ~astropy.units.quantity.Quantity

        """
        self.r = r
        self.theta = theta
        self.phi = phi
        self.a = a
        super().__init__(r.si.value, theta.si.value, phi.si.value, a=a.si.value)
        self.system = "BoyerLindquist"
        self._dimension = {
            "r": self.r,
            "theta": self.theta,
            "phi": self.phi,
            "a": self.a,
            "system": self.system,
        }
        self._dimension_order = ("r", "theta", "phi")

    def __repr__(self):
        return "Boyer-Lindquist r: {}, theta: {}, phi: {} | a: {}".format(
            self.r, self.theta, self.phi, self.a
        )

    def __str__(self):
        return self.__repr__()

    def __getitem__(self, item):
        """
        Method to return coordinates.
        Objects are subsctiptable with both explicit names of parameters
        and integer indices.

        Parameters
        ----------
        item : str or int
                Name of the parameter or its index.
                If ``'system'``, Name of coordinate is returned.
                If ``'a'``, spin factor of the body, ``self.a`` is returned.
        """
        if isinstance(item, (int, np.integer)):
            return self._dimension[self._dimension_order[item]]
        return self._dimension[item]

    def si_values(self):
        """
        Function for returning values in SI units.

        Returns
        -------
        ~numpy.ndarray
            Array containing values in SI units (m, rad, rad)

        """
        element_list = [self.r.to(u.m), self.theta.to(u.rad), self.phi.to(u.rad)]
        return np.array([e.value for e in element_list], dtype=float)

    def to_cartesian(self):
        """
        Method for conversion to cartesian coordinates.

        Returns
        -------
        ~einsteinpy.coordinates.core.Cartesian
            Cartesian representation of the BL Coordinates.

        """
        x, y, z = self.convert_cartesian()
        return Cartesian(x * u.m, y * u.m, z * u.m)

    def to_spherical(self):
        """
        Method for conversion to spherical coordinates.

        Returns
        -------
        ~einsteinpy.coordinates.core.Spherical
            Spherical representation of the BL Coordinates.

        """
        r, theta, phi = self.convert_spherical()
        return Spherical(r * u.m, theta * u.rad, phi * u.rad)
