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
        self.x_u = u.m.to(x.unit, 1)
        self.y_u = u.m.to(y.unit, 1)
        self.z_u = u.m.to(z.unit, 1)
        self.x = x.si.value
        self.y = y.si.value
        self.z = z.si.value
        self.system = "Cartesian"

    def __repr__(self):
        return "Cartesian x: {}, y: {}, z: {}".format(self.x * self.x_u, self.y * self.y_u, self.z * self.z_u)

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
        element_list = [self.x, self.y, self.z]
        return np.array([e for e in element_list], dtype=float)

    def norm(self):
        """
        Function for finding euclidean norm of a vector.

        Returns
        -------
        ~astropy.units.quantity.Quantity
            Euclidean norm with units.

        """
        return (np.sqrt(self.x ** 2 + self.y ** 2 + self.z ** 2) * u.m)

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

        return (x + y + z) * u.m

    def to_spherical(self):
        """
        Method for conversion to spherical coordinates.

        Returns
        -------
        ~einsteinpy.coordinates.core.Spherical
            Spherical representation of the Cartesian Coordinates.

        """
        r = self.norm().value
        theta = np.arccos(self.z / r)
        phi = np.arctan2(self.y, self.x)

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
        a = a.to(u.m).value
        w = self.norm().value ** 2 - a ** 2
        r = np.sqrt(0.5 * (w + np.sqrt((w ** 2) + (4 * (a ** 2) * (self.z ** 2)))))
        theta = np.arccos(self.z / r)
        phi = np.arctan2(self.y, self.x)

        return BoyerLindquist(r * u.m, theta * u.rad, phi * u.rad, a * u.m)


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
        self.r_u = u.m.to(r.unit, 1)
        self.theta_u = u.rad.to(theta.unit, 1)
        self.phi_u = u.rad.to(phi.unit, 1)
        self.r = r.si.value
        self.theta = theta.si.value
        self.phi = phi.si.value
        self.system = "Spherical"

    def __repr__(self):
        return "Spherical r: {}, theta: {}, phi: {}".format(
            self.r * self.r_u, self.theta * self.theta_u, self.phi * self.phi_u
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
        element_list = [self.r, self.theta, self.phi]
        return np.array([e for e in element_list], dtype=float)

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
        self.r_u = u.m.to(r.unit, 1)
        self.theta_u = u.rad.to(theta.unit, 1)
        self.phi_u = u.rad.to(phi.unit, 1)
        self.a_u = u.m.to(a.unit, 1)
        self.r = r.si.value
        self.theta = theta.si.value
        self.phi = phi.si.value
        self.a = a.si.value
        self.system = "BoyerLindquist"

    def __repr__(self):
        return "Boyer-Lindquist r: {}, theta: {}, phi: {} | a: {}".format(
            self.r * self.r_u, self.theta * self.theta_u, self.phi * self.phi_u, self.a * self.a_u
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
        element_list = [self.r, self.theta, self.phi]
        return np.array([e for e in element_list], dtype=float)

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

        return Cartesian(x * u.m, y * u.m, z * u.m)

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
