import astropy.units as u
import numpy as np

# utils may not be required
from einsteinpy import utils # - ????


from einsteinpy import constant, utils
from einsteinpy.coordinates import BoyerLindquist, KerrSchild, Spherical
# from einsteinpy.utils import schwarzschild_radius_dimensionless


_c = constant.c
_G = constant.G
_Cc = constant.coulombs_const

class Metric:
    """
    Class for defining general Metric Tensors

    A place to gather functions from einsteinpy.utils - ?????

    This module + `utils.scalar_factor` (with perhaps a name change) 
    should have all relevant Numerical Relativity utilities

    Main usage will be in perturbative treatment 
    of EFE and its solutions

    @staticmethods: to maintain their usage as utility functions,
    even when they are logically connected to this class - ?????
    """

    # Precomputed list of tuples consisting of indices 
    # of christoffel symbols which are non-zero in the Metric
    # Purpose: To exploit symmetries to reduce number of computations
    nonzero_christoffels_list = []

    def __init__(self, arr, coords, M, a, Q, name="GenericMetricTensor"):
        """
        Constructor

        Units: Geometrized Unit System
        
        Parameters
        ----------
        arr : ~numpy.array
            Numpy Array containing Metric components
        coords : einsteinpy.coordinates
            Coordinate system, in which Metric has been represented
        M : float
            Mass of gravitating body, e.g. Black Hole
        a : float
            Spin Parameter, 0 <= a <= 1 - ???? (Only if Geom Units)
        Q : float
            Charge on gravitating body, e.g. Black Hole
        name : str
            Name of the Metric Tensor. Defaults to "GenericMetricTensor"
        """
        self.name = name
        # - ?????
        # Need to decide how to allow arbitray metric,
        # that connects well with `coords`,
        # and leaves predefined metric classes immutable
        self.arr = arr
        self.coords = coords
        self.M = M
        self.a = a
        self.Q = Q


    def __str__(self): # - ????? (arr)
        return (
            " ( name: ({0}), coordinates: ({1}), array: ({2}). mass: ({3}), spin parameter: ({4}), charge: ({"
            "5}) )".format(
                self.name, self.coords, self.arr, self.M, self.a, self.Q
            )
        )

    def __repr__(self): # - ????? (arr)
        return (
            " ( name: ({0}), coordinates: ({1}), array: ({2}). mass: ({3}), spin parameter: ({4}), charge: ({"
            "5}) )".format(
                self.name, self.coords, self.arr, self.M, self.a, self.Q
            )
        )

    @staticmethod
    def sigma(r, theta, a):
        """
        Returns the value r^2 + a^2 * cos^2(theta)
        Specific to Boyer-Lindquist coordinates

        Parameters
        ----------
        r : float
            r-component of 4-Position
        theta : float
            theta-component of 4-Position
        a : float
            Spin Parameter

        Returns
        -------
        float
            The value of ``r^2 + a^2 * cos^2(theta)``
        """
        sigma = (r ** 2) + ((a * np.cos(theta)) ** 2)
        return sigma

    @staticmethod
    def delta(r, M, a, Q=0):
        """
        Returns the value of (r^2 - Rs * r + a^2 + Rq^2)
        Specific to Boyer-Lindquist coordinates
        Applies to Kerr Geometry

        Parameters
        ----------
        r : float
            r-component of 4-Position
        M : float
            Mass of gravitating body
        a : float
            Spin Parameter
        Q : float
            Charge on gravitating body
            Defaults to 0 (for Kerr Geometry)

        Returns
        -------
        float
            The value of ``r^2 - Rs * r + a^2 + Rq^2``

        """
        # Removing c, G, Cc as parameters - ?????
        Rs = 2 * M * _G.value / _c.value ** 2
        # Square of Geometrized Charge
        RQ2 = (Q ** 2) * _G.value * _Cc.value / _c.value ** 4
        delta = (r ** 2) - (Rs * r) + (a ** 2) + RQ2
        return delta

    @staticmethod
    def rho(r, theta, a):
        """
        Returns the value of sqrt(r^2 + a^2 * cos^2(theta)) == sqrt(sigma)
        Specific to Boyer-Lindquist coordinates
        Applies to Kerr-Newman Geometry

        Parameters
        ----------
        r : float
            r-component of 4-Position
        theta : float
            theta-component of 4-Position
        a : float
            Spin Parameter

        Returns
        -------
        float
            The value of ``sqrt(r^2 + a^2 * cos^2(theta))``

        """
        return np.sqrt(Metric.sigma(r, theta, a))

    @staticmethod
    def r_ks(x, y, z, a):
        """
        Returns the value of r, after solving (x**2 + y**2) / (r**2 + a**2) + z**2 / r**2 = 1
        'r' is not the Radius Coordinate of Spherical Polar or BL Coordinates
        Specific to Cartesian form of Kerr-Schild Coordinates
        """
        pass

    @staticmethod
    def schwarzschild_radius(M):
        """
        Returns Schwarzschild Radius in SI units

        Parameters
        ----------
        M : ~astropy.units.kg
            Mass of gravitating body, in kg
        Returns
        -------
        r : ~astropy.units.m
            Schwarzschild Radius for a given mass, in m
        """
        if not isinstance(M, u.quantity.Quantity):
            M = M * u.kg
        # Removing c & G as parameters - ?????
        M = M.to(u.kg)
        num = 2 * _G * M
        deno = _c ** 2
        return num / deno

    @staticmethod
    def schwarzschild_radius_dimensionless(M):
        """
        Returns the value of Schwarzschild Radius

        Parameters
        ----------
        M : float
            Mass of gravitating body
        Returns
        -------
        Rs : float
            Schwarzschild Radius for a given mass
        """
        # Removing c & G as parameters - ?????
        # Unused - ?????
        Rs = 2 * M * _G.value / _c.value ** 2
        return Rs

    @staticmethod
    def spin_parameter(J, M):
        """
        Returns Spin Parameter (a) of a Rotating Body, in SI units

        Parameters
        ----------
        J : float
            Angular momentum, in SI units
        M : float
            Mass of body, in SI units

        Returns
        -------
        float
            Spin Parameter

        """
        # Removing c as parameter - ?????
        return J / (M * _c)
    
    @staticmethod
    def scaled_spin_parameter(a, M):
        """
        Returns scaled Spin Parameter (a), to incorporate changed units
        
        Parameters
        ----------
        a : float
            Number between 0 and 1
        M : float
            Mass of gravitating body
        
        Returns
        -------
        float
            Scaled Spin Parameter

        Raises
        ------
        ValueError
            If a is not between 0 and 1

        """
        # Where is this to be used - ?????
        # What are we scaling to - ????? (Mass?)
        # Can perhaps put in einsteinpy.units
        # Removing c & G as parameters - ?????
        half_rs = _G.value * M / _c.value ** 2
        if a < 0 or a > 1:
            raise ValueError("a should be between 0 and 1.")
        return a * half_rs

    @staticmethod
    def charge_geometrized(Q):
        """
        Geometrized representation of the Electric Charge on the gravitating body

        Parameters
        ----------
        Q : float
            Charge on gravitating body

        Returns
        -------
        float
            Geometrized Charge

        """
        # Removing c, G, Cc as parameters - ?????
        # Belongs in the `einsteinpy.units` module - ?????
        # Unused - ?????
        return (Q / (_c.value ** 2)) * np.sqrt(_G.value * _Cc.value)
    
    @staticmethod
    def singularities(M, a, Q=0, coords="BL"):
        """
        Returns the Singularities of the Metric
        Depends on the choice of Coordinate Systems
        Applies to Kerr and Kerr-Newman Geometries

        Parameters
        ----------
        M : float
            Mass of gravitating body
        a : float
            Spin Parameter
        Q : float
            Charge on gravitating body in the Metric
            Defaults to 0 (for Kerr Geometry)
        coords : str
            Coordinates, in which to calculate and return the singularities
            "BL" for Boyer-Lindquist, "KS" for Kerr-Schild
            Defaults to "BL"

        Returns
        -------
        dict - ????? (Docs)
            Dictionary of singularities in the geometry
            ``{
                "inner_ergosphere": lambda function
                "inner_horizon": float
                "outer_horizon": float
                "outer_ergosphere": lambda function
            }``

        """
        # Removing c, G, Cc as parameters - ?????
        Rs = 2 * M * _G.value / _c.value ** 2
        # Square of Geometrized Charge
        RQ2 = (Q ** 2) * _G.value * _Cc.value / _c.value ** 4

        inner_ergosphere = None
        inner_horizon = 0
        outer_horizon = 0
        outer_ergosphere = None
        
        if coords == "BL":
            inner_ergosphere = lambda theta: (Rs - np.sqrt((Rs ** 2) - (4 * (a * np.cos(theta)) ** 2) - (4 * RQ2))) / 2
            inner_horizon = (Rs - np.sqrt((Rs ** 2) - (4 * a ** 2) - (4 * RQ2))) / 2
            outer_horizon = (Rs + np.sqrt((Rs ** 2) - (4 * a ** 2) - (4 * RQ2))) / 2
            outer_ergosphere = lambda theta: (Rs + np.sqrt((Rs ** 2) - (4 * (a * np.cos(theta)) ** 2) - (4 * RQ2))) / 2
        
        elif coords == "KS":
            # - ????? (To be filled in, after refactoring `coordinates`)
            pass

        return {
                "inner_ergosphere": inner_ergosphere,
                "inner_horizon": inner_horizon,
                "outer_horizon": outer_horizon,
                "outer_ergosphere": outer_ergosphere,
            }


    # The methods below are based on how a user wants to set up a problem
    # Potentially, this should lead to support for exotic geometries
    # Some functions below are to be used as blueprints for
    # defining classes for well-known Metric Tensors.
    def metric(self, r, theta, M, a):
        """
        Returns the Metric in its general form
        """
        pass

    def metric_inv(self, r, theta, M, a):
        """
        Returns the inverse of Metric
        """
        pass

    def dmetric_dx(self, r, theta, M, a):
        """
        Returns differentiation of each Metric component w.r.t. coordinates
        """
        dmdx = np.zeros(shape=(4, 4, 4), dtype=float)
        def due_to_x0():
            # Diff wrt x_0 component
            nonlocal dmdx
            pass

        def due_to_x1():
            # Diff wrt x_1 component
            nonlocal dmdx
            pass

        def due_to_x2():
            # Diff wrt x_2 component
            nonlocal dmdx
            pass

        def due_to_x3():
            # Diff wrt x_3 component
            nonlocal dmdx
            pass

        due_to_x0()
        due_to_x1()
        due_to_x2()
        due_to_x3()
        return dmdx

    def christoffels(self, r, theta, M, a):
        """
        Returns Christoffel Symbols for the Metric
        """
        pass

    def nonzero_christoffels(self):
        """
        Returns a list of tuples consisting of indices 
        of christoffel symbols, which are non-zero
        in the Metric, computed in real-time.
        """
        pass
    """ FUNCTIONS TO INCLUDE -- (END) - ???? """




    """ FUNCTIONS THAT MAY BE INCLUDED HERE - ???? """
    def em_potential_covariant(self, r, theta, M, a, Q):
        """
        Returns Electromagnetic 4-Potential
        Specific to Kerr-Newman Geometry

        Check Eq. (1), Page 2: https://arxiv.org/pdf/1407.1530.pdf - ?????
        Or, Wikipedia: https://en.wikipedia.org/wiki/Kerr%E2%80%93Newman_metric#Electromagnetic_field_tensor_in_Boyer-Lindquist_form

        Parameters
        ----------
        r : float
            r-component of 4-Position
        theta : float
            theta-component of 4-Position
        M : float
            Mass of gravitating body
        a : float
            Spin Parameter
        Q : float
            Charge on gravitating body
        
        Returns
        -------
        ~numpy.array
            Electromagnetic 4-Potential
            Numpy array of shape (4,)

        """
        # Removing c, G, Cc as parameters - ?????
        # Square of Geometrized Charge
        RQ = np.sqrt((Q ** 2 * _G * _Cc) / _c ** 4)
        rho2 = Metric.rho(r, theta, a) ** 2
        pot = np.zeros((4,), dtype=float)
        pot[0] = r * RQ / rho2
        # (_c ** 2 / _G * M) is extraneous - ?????
        pot[3] = (_c ** 2 / _G * M) *  (-r * a * RQ * np.sin(theta) ** 2 / rho2)

        return pot

    def maxwell_tensor_covariant(self, r, theta, M, a, Q):
        """
        Returns Covariant Maxwell Stress Tensor
        Specific to Kerr-Newman Geometry

        Parameters
        ----------
        r : float
            r-component of 4-Position
        theta : float
            theta-component of 4-Position
        M : float
            Mass of gravitating body
        a : float
            Spin Parameter
        Q : float
            Charge on gravitating body
        
        Returns
        -------
        ~numpy.array
            Maxwell Stress Tensor
            Numpy array of shape (4, 4)

        """
        # Removing c, G, Cc as parameters - ?????
        pass

    def maxwell_tensor_contravariant( r, theta, M, a, Q):
        """
        Returns Contravariant Maxwell Stress Tensor
        Specific to Kerr-Newman Geometry
        """
        pass
    """ FUNCTIONS THAT MAY BE INCLUDED HERE -- (END) - ???? """


    """ NOT SURE, IF THESE SHOULD BE HERE - ???? """
    # @u.quantity_input(mass=u.kg)
    def time_velocity(pos_vec, vel_vec, mass, a):
        """
        # Velocity of coordinate time wrt proper metric
        Timelike component of 4-Velocity

        Parameters
        ----------
        pos_vector : ~numpy.array
            Vector with r, theta, phi components in SI units
        vel_vector : ~numpy.array
            Vector with velocities of r, theta, phi components in SI units
        mass : ~astropy.units.kg
            Mass of gravitating body
        a : float
            Any constant
        Q : ~astropy.units.C
        Charge on gravitating body

        Returns
        -------
        ~astropy.units.one
            # Velocity of time
            Timelike component of 4-Velocity

        """
        # Similar function for Kerr and Kerr-Newman
        # Perhaps, this should be in coordinates in some form
        pass
    """ NOT SURE, IF THESE SHOULD BE HERE -- (END) - ???? """
