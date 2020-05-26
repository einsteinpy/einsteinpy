import numpy as np

# utils may not be required
from einsteinpy import utils # - ????
from einsteinpy.constant import * # - ????
from einsteinpy.coordinates import * # - ????


class Metric:
    """
    Class for defining general Metric Tensors

    A place to gather functions from einsteinpy.utils - ????

    This module + `utils.scalar_factor` (with perhaps a name change) 
    should have all relevant Numerical Relativityutilities

    Main usage will be in perturbative treatment 
    of EFE and its solutions
    """

    # Precomputed list of tuples consisting of indices 
    # of christoffel symbols which are non-zero in the Metric
    # Purpose: To exploit symmetries to reduce number of computations
    nonzero_christoffels_list = []

    def __init__(self, arr, coords, M, a, q, name):
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
            Mass of the gravitating object, e.g. Black Hole
        a : float
            Spin Parameter 0 <= a <= 1
        q : float
            Charge
        name : str
            Name of the Metric Tensor. Defaults to "GenericMetricTensor".
        """
        self.name = name
        self.arr = arr
        self.coords = coords
        self.M = M
        self.a = a
        self.q = q


    def __str__(self):
        pass

    def __repr__(self):
        pass


    """ FUNCTIONS TO NOT INCLUDE, WITH REASONS - ???? """
    def scaled_spin_factor(a, M, c=constant.c.value, G=constant.G.value):
        """
        Returns a scaled version of spin factor(a) ❌
        """
        # Pivoting to Geom
        pass

    def schwarzschild_radius(mass, c=constant.c, G=constant.G):
        """
        Returns Schwarzschild Radius (In SI) ❌
        """
        # Pivoting to Geom
        # Combining two functions into 1
        pass

    def schwarzschild_radius_dimensionless(M, c=constant.c.value, G=constant.G.value):
        """
        Returns Schwarzschild Radius (Dimensionless) ❌
        """
        # Pivoting to Geom
        # Combining two functions into 1
        pass
    """ FUNCTIONS TO NOT INCLUDE, WITH REASONS -- (END) - ???? """


    # FUNCTIONS TO INCLUDE - ????
    @staticmethod
    def sigma(r, theta, a):
        """
        Returns the value r^2 + a^2 * cos^2(theta)
        Specific to Boyer-Lindquist coordinates
        """
        # Specific to BL, but can be included as a @staticmethod 
        pass

    @staticmethod
    def delta(r, M, a):
        """
        Returns the value r^2 - Rs * r + a^2
        Specific to Boyer-Lindquist coordinates
        """
        # Specific to BL, but can be included as a @staticmethod 
        pass
    
    @staticmethod
    def r_ks(x, y, z, a):
        """
        Returns the value of r, after solving (x**2 + y**2) / (r**2 + a**2) + z**2 / r**2 = 1
        'r' is not the Radius Coordinate of Spherical Polar or BL Coordinates
        Specific to Cartesian form of Kerr-Schild Coordinates
        """
        pass

    @staticmethod
    def rho(r, theta, a):
        """
        Returns the value sqrt(r^2 + a^2 * cos^2(theta)).
        Specific to Boyer-Lindquist coordinates
        Specific to Kerr-Newman Geometry
        """
        pass

    @staticmethod
    def delta(r, M, a, Q):
        """
        Returns the value r^2 - Rs * r + a^2
        Specific to Boyer-Lindquist coordinates
        Specific to Kerr-Newman Geometry
        """
        pass


    @staticmethod
    def spin_factor(J, M):
        """
        Calculate spin factor(a) of kerr body
        """
        pass

    @staticmethod
    def schwarzschild_radius(M): # _dimensionless
        """
        Returns Schwarzschild Radius
        """
        pass

    @staticmethod
    def event_horizon(M, a, theta=np.pi / 2, coord="BL"):
        """
        Calculate the radius of event horizon
        """
        pass

    @staticmethod
    def radius_ergosphere(M, a, theta=np.pi / 2, coord="BL",):
        """
        Calculate the radius of ergospere ~of Kerr black hole~ at a specific azimuthal angle
        """
        pass

    @staticmethod
    def charge_length_scale(Q):
        """
        Returns a length scale corrosponding to the Electric Charge Q of the mass
        Specific to Kerr-Newman Geometry
        """
        pass


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
    def em_potential(r, theta, a, Q, M):
        """
        Returns 4-Potential
        Specific to Kerr-Newman Geometry
        """
        pass

    def maxwell_tensor_covariant( r, theta, a, Q, M):
        """
        Returns Covariant Maxwell Stress Tensor
        """
        pass

    def maxwell_tensor_contravariant( r, theta, a, Q, M):
        """
        Returns Contravariant Maxwell Stress Tensor
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
            Mass of the body
        a : float
            Any constant
        Q : ~astropy.units.C
        Charge on the massive body

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
