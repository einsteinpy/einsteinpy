import astropy.units as u
import numpy as np


class Cartesian:
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
        self.system = "Cartesian"

    def __repr__(self):
        return "Cartesian x: {}, y: {}, z: {}".format(self.x, self.y, self.z)

    def __str__(self):
        return self.__repr__()

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
        r = self.norm()
        theta = np.arccos(self.z / r)
        phi = np.arctan2(self.y, self.x)

        return Spherical(r, theta, phi)

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
        w = self.norm() ** 2 - a ** 2
        r = np.sqrt(0.5 * (w + np.sqrt((w ** 2) + (4 * (a ** 2) * (self.z ** 2)))))
        theta = np.arccos(self.z / r)
        phi = np.arctan2(self.y, self.x)

        return BoyerLindquist(r, theta, phi, a)


class Spherical:
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
        self.system = "Spherical"

    def __repr__(self):
        return "Spherical r: {}, theta: {}, phi: {}".format(
            self.r, self.theta, self.phi
        )

    def __str__(self):
        return self.__repr__()

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
        x = self.r * np.cos(self.phi) * np.sin(self.theta)
        y = self.r * np.sin(self.phi) * np.sin(self.theta)
        z = self.r * np.cos(self.theta)

        return Cartesian(x, y, z)

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
        cart = self.to_cartesian()
        return cart.to_bl(a)


class BoyerLindquist:
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
        self.system = "BoyerLindquist"

    def __repr__(self):
        return "Boyer-Lindquist r: {}, theta: {}, phi: {} | a: {}".format(
            self.r, self.theta, self.phi, self.a
        )

    def __str__(self):
        return self.__repr__()

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
        sin_norm = np.sqrt(self.r ** 2 + self.a ** 2) * np.sin(self.theta)
        x = sin_norm * np.cos(self.phi)
        y = sin_norm * np.sin(self.phi)
        z = self.r * np.cos(self.theta)

        return Cartesian(x, y, z)

    def to_spherical(self):
        """
        Method for conversion to spherical coordinates.

        Returns
        -------
        ~einsteinpy.coordinates.core.Spherical
            Spherical representation of the BL Coordinates.

        """
        cart = self.to_cartesian()
        return cart.to_spherical()
