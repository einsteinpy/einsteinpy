import astropy.units as u
import numpy as np

from einsteinpy import constant

_c = constant.c
_G = constant.G
_Cc = constant.coulombs_const


class Metric:
    """
    Class for defining general Metric Tensors

    To be used in perturbative treatment of EFE and its solutions
   
    """

    @u.quantity_input(M=u.kg, a=u.km, Q=u.C)
    def __init__(
        self,
        coords,
        M,
        a=0,
        Q=0,
        name="Generic Metric",
        metric_cov=None,
        christoffels=None,
        f_vec=None,
        perturbation=None,
    ):
        """
        Constructor
        
        Parameters
        ----------
        coords : string
            Coordinate system, in which Metric is to be represented
            "S" - Schwarzschild: Only applicable to Schwarzschild solutions
            "BL" - Boyer-Lindquist: Applicable to Kerr-Newman solutions
            "KS" - Kerr-Schild: Useful for adding perturbations to Kerr-Newman solutions
        M : float
            Mass of gravitating body, e.g. Black Hole
        a : float
            Spin Parameter
            Defaults to ``0``
        Q : float
            Charge on gravitating body, e.g. Black Hole
            Defaults to ``0``
        name : str
            Name of the Metric Tensor. Defaults to ``"Generic Metric"``
        metric_cov : function
            Function, defining Covariant Metric Tensor
            It should return a real-valued tensor (2D Array), at supplied coordinates
            Defaults to ``None``
            Consult pre-defined classes for function definition
        christoffels : function
            Function, defining Christoffel Symbols
            It should return a real-valued (4,4,4) array, at supplied coordinates
            Defaults to ``None``
            Consult pre-defined classes for function definition
        f_vec : function
            Function, defining RHS of Geodesic Equation
            It should return a real-valued (8) vector, at supplied coordinates
            Defaults to ``None``
            Consult pre-defined classes for function definition
        perturbation : function
            Function, defining Covariant Perturbation Tensor
            It should return a real-valued tensor (2D Array), at supplied coordinates
            Defaults to ``None``
            Function definition similar to ``metric_cov``
        
        """
        self.name = name
        self.coords = coords
        self.M = M
        self.a = a
        self.Q = Q
        self.metric_cov = metric_cov
        self.christoffels = christoffels
        self.f_vec = f_vec
        self.perturbation = perturbation
        # Need physical checks for `perturbation`
        # Expert opinion needed

        self.sch_rad = self.schwarzschild_radius(M)

    def __str__(self):
        return f"( Name: ({self.name}), Coordinates: ({self.coords}), \
            Mass: ({self.M}), Spin parameter: ({self.a}), \
            Charge: ({self.Q}), Schwarzschild Radius: ({self.sch_rad}) )"

    def __repr__(self):
        return f"( Name: ({self.name}), Coordinates: ({self.coords}), \
            Mass: ({self.M}), Spin parameter: ({self.a}), \
            Charge: ({self.Q}), Schwarzschild Radius: ({self.sch_rad}) )"

    @staticmethod
    @u.quantity_input(r=u.km, theta=u.rad, a=u.km)
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
    @u.quantity_input(r=u.km, M=u.kg, a=u.km, Q=u.C)
    def delta(r, M, a, Q=0):
        """
        Returns the value of (r^2 - r_s * r + a^2 + r_Q^2)
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
            Defaults to ``0`` (for Kerr Geometry)

        Returns
        -------
        float
            The value of ``r^2 - r_s * r + a^2 + r_Q^2``

        """
        r_s = 2 * M * _G.value / _c.value ** 2
        # Square of Geometrized Charge
        r_Q2 = (Q ** 2) * _G.value * _Cc.value / _c.value ** 4
        delta = (r ** 2) - (r_s * r) + (a ** 2) + r_Q2
        return delta

    @staticmethod
    @u.quantity_input(r=u.km, theta=u.rad, a=u.km)
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
    @u.quantity_input(x=u.km, y=u.km, z=u.km, a=u.km)
    def r_ks(x, y, z, a):
        """
        Returns the value of r, after solving (x**2 + y**2) / (r**2 + a**2) + z**2 / r**2 = 1
        'r' is not the Radius Coordinate of Spherical Polar or BL Coordinates
        Specific to Cartesian form of Kerr-Schild Coordinates
        """
        raise NotImplementedError

    @staticmethod
    @u.quantity_input(M=u.kg)
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
        return 2 * _G * M / _c ** 2

    @staticmethod
    @u.quantity_input(J=u.kg * u.m ** 2 / u.s, M=u.kg)
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
        return J / (M * _c)

    @staticmethod
    @u.quantity_input(a=u.km, M=u.kg)
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
        # Unused - Not used in any method ?????
        # What are we scaling to - ????? (Mass?)
        # Can perhaps put in einsteinpy.units?
        half_rs = _G.value * M / _c.value ** 2
        if a < 0 or a > 1:
            raise ValueError("a should be between 0 and 1.")
        return a * half_rs

    @staticmethod
    @u.quantity_input(Q=u.C)
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
        # Unused - Not used in any method ?????
        # What are we scaling to - ????? (Mass?)
        # Can perhaps put in einsteinpy.units?
        return (Q / (_c.value ** 2)) * np.sqrt(_G.value * _Cc.value)

    @staticmethod
    @u.quantity_input(M=u.kg, a=u.km, Q=u.C)
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
            Defaults to ``0 `(for Kerr Geometry)
        coords : str
            Coordinate System, in which singularities are calculated
            "BL" for Boyer-Lindquist, "KS" for Kerr-Schild
            Defaults to ``"BL"``

        Returns
        -------
        dict
            Dictionary of singularities in the geometry
            ``{
            "inner_ergosphere": function(theta),
            "inner_horizon": float,
            "outer_horizon": float,
            "outer_ergosphere": function(theta)
            }``

        """
        r_s = 2 * M * _G.value / _c.value ** 2
        # Square of Geometrized Charge
        r_Q2 = (Q ** 2) * _G.value * _Cc.value / _c.value ** 4

        def _in_ergo(theta):
            return (
                r_s
                - np.sqrt((r_s ** 2) - (4 * (a * np.cos(theta)) ** 2) - (4 * r_Q2)) / 2
            )

        def _out_ergo(theta):
            return (
                r_s
                + np.sqrt((r_s ** 2) - (4 * (a * np.cos(theta)) ** 2) - (4 * r_Q2)) / 2
            )

        if coords == "S":  # Schwarzschild Geometry
            return {
                "inner_ergosphere": 0,
                "inner_horizon": 0,
                "outer_horizon": r_s,
                "outer_ergosphere": r_s,
            }

        elif coords == "BL":  # Kerr & Kerr-Newman Geometries
            return {
                "inner_ergosphere": _in_ergo,
                "inner_horizon": (r_s - np.sqrt((r_s ** 2) - (4 * a ** 2) - (4 * r_Q2)))
                / 2,
                "outer_horizon": (r_s + np.sqrt((r_s ** 2) - (4 * a ** 2) - (4 * r_Q2)))
                / 2,
                "outer_ergosphere": _out_ergo,
            }

        elif coords == "KS":  # Kerr & Kerr-Newman Geometries
            # - ????? (To be filled in, after refactoring `coordinates`)
            raise NotImplementedError

    # Derived classes should only define metric_covariant() function
    # Check Kerr or Kerr Newman for understanding this
    # No need for metric_contravariant()
    def metric_covariant(self, x_vec):
        """
        Returns Covariant Metric Tensor
        Adds Kerr-Schild (Linear) Perturbation to metric, \
        if ``perturbation`` is defined in Metric object

        Parameters
        ----------
        x_vec : numpy.array
            Position 4-Vector

        Returns
        -------
        ~numpy.array
            Covariant Metric Tensor
            Numpy array of shape (4,4)

        """
        g_cov = self.metric_cov(x_vec, self.M, self.a, self.Q)

        if self.perturbation:
            p_cov = self.perturbation(x_vec, self.M, self.a, self.Q)
            return g_cov + p_cov

        return g_cov

    def metric_contravariant(self, x_vec):
        """
        Returns Contravariant Metric Tensor
        Adds Kerr-Schild (Linear) Perturbation to metric, \
        if ``perturbation`` is not None in Metric object
        
        Parameters
        ----------
        x_vec : numpy.array
            Position 4-Vector
        
        Returns
        -------
        ~numpy.array
            Contravariant Metric Tensor

        """
        g_cov = self.metric_covariant(x_vec)
        g_contra = np.linalg.inv(g_cov)

        return g_contra

    # Christoffels to be supplied as a function
    # f_vec to be supplied as a function

    # M not needed as parameter in 5 functions below - ????? \
    # (if (_c.value ** 2 / _G.value * M) is erroneous)
    # Check Eq. (1), Page 2: https://arxiv.org/pdf/1407.1530.pdf - ?????
    # Or: https://en.wikipedia.org/wiki/Kerr%E2%80%93Newman_metric#Electromagnetic_field_tensor_in_Boyer-Lindquist_form
    @u.quantity_input(r=u.km, theta=u.rad, M=u.kg, a=u.km, Q=u.C)
    def em_potential_covariant(self, r, theta, M, a, Q):
        """
        Returns Covariant Electromagnetic 4-Potential
        Specific to Kerr-Newman Geometries

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
            Covariant Electromagnetic 4-Potential
            Numpy array of shape (4,)

        """
        # Square of Geometrized Charge
        r_Q = np.sqrt((Q ** 2 * _G.value * _Cc.value) / _c.value ** 4)
        rho2 = Metric.rho(r, theta, a) ** 2

        A = np.zeros((4,), dtype=float)
        A[0] = r * r_Q / rho2
        # (_c ** 2 / _G * M) is extraneous - ?????
        A[3] = (_c.value ** 2 / _G.value * M) * (
            -r * a * r_Q * np.sin(theta) ** 2 / rho2
        )

        return A

    @u.quantity_input(r=u.km, theta=u.rad, M=u.kg, a=u.km, Q=u.C)
    def em_potential_contravariant(self, r, theta, M, a, Q):
        """
        Returns Contravariant Electromagnetic 4-Potential
        Specific to Kerr-Newman Geometries

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
            Contravariant Electromagnetic 4-Potential
            Numpy array of shape (4,)

        """
        A_cov = self.em_potential_covariant(r, theta, M, a, Q)
        x_vec = [0, r, theta, 0]  # t & phi have no bearing on Metric
        g_contra = self.metric_contravariant(x_vec=x_vec)
        # @ has similar perf to np.dot or np.matmul
        # https://stackoverflow.com/a/58116209/11922029
        return g_contra @ A_cov

    @u.quantity_input(r=u.km, theta=u.rad, M=u.kg, a=u.km, Q=u.C)
    def maxwell_tensor_covariant(self, r, theta, M, a, Q):
        """
        Returns Covariant Maxwell Stress Tensor
        Specific to Kerr-Newman Geometries

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
            Covariant Maxwell Stress Tensor
            Numpy array of shape (4, 4)

        """
        r_Q = np.sqrt((Q ** 2 * _G.value * _Cc.value) / _c.value ** 4)
        rho2 = Metric.rho(r, theta, a) ** 2
        # Partial derivatives of rho2
        drho2_dr = 2 * r
        drho2_dtheta = -(a ** 2 * np.sin(2 * theta))

        F = np.zeros((4, 4), dtype=float)

        F[0, 1] = -(r_Q * (rho2 - drho2_dr * r)) / (rho2 ** 2)
        F[1, 0] = -F[0, 1]
        F[0, 2] = (r * r_Q * drho2_dtheta) / (rho2 ** 2)
        F[2, 0] = -F[0, 2]
        # (_c.value ** 2 / _G.value * M) is extraneous - ????? (Issue #144, perhaps)
        F[1, 3] = (
            (_c.value ** 2 / _G.value * M)
            * (1 / rho2 ** 2)
            * (a * r_Q * np.sin(theta) ** 2)
            * (rho2 ** 2 - 2 * r ** 2)
        )
        F[3, 1] = -F[1, 3]
        # (_c.value ** 2 / _G.value * M) is extraneous - ????? (Issue #144, perhaps)
        F[2, 3] = (
            (_c.value ** 2 / _G.value * M)
            * (1 / rho2 ** 2)
            * (a * r_Q * r * np.sin(2 * theta))
            * (rho2 + (a * np.sin(theta) ** 2))
        )
        F[3, 2] = -F[2, 3]

        return F

    @u.quantity_input(r=u.km, theta=u.rad, M=u.kg, a=u.km, Q=u.C)
    def maxwell_tensor_contravariant(self, r, theta, M, a, Q):
        """
        Returns Contravariant Maxwell Stress Tensor
        Specific to Kerr-Newman Geometries

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
            Contravariant Maxwell Stress Tensor
            Numpy array of shape (4, 4)

        """
        F_cov = self.maxwell_tensor_covariant(r, theta, M, a, Q)
        x_vec = [0, r, theta, 0]
        g_contra = self.metric_contravariant(x_vec=x_vec)

        F_contra = g_contra @ F_cov @ g_contra

        return F_contra

    # Functions below are untouched
    @staticmethod
    @u.quantity_input(t=u.s)
    def scalar_factor(t, era="md", tuning_param=1.0):
        """
        Acceleration of the universe in cosmological models of Robertson Walker
        Flat Universe.

        Parameters
        ----------
        era : string
            Can be chosen from 'md' (Matter Dominant),
            'rd' (Radiation Dominant) and 'ded' (Dark Energy Dominant)
        t : ~astropy.units.s
            Time for the event
        tuning_param : float, optional
            Unit scaling factor, defaults to 1

        Returns
        -------
        float
            Value of scalar factor at time t.

        Raises
        ------
        ValueError : If era is not 'md' , 'rd', and 'ded'.

        """
        T = t.to(u.s).value
        if era == "md":
            return tuning_param * (T ** (2 / 3))
        elif era == "rd":
            return tuning_param * (T ** (0.5))
        elif era == "ded":
            hubble_const = (constant.Cosmo_Const / 3) ** 0.5
            val = np.e ** (hubble_const.value * T)
            return tuning_param * val
        else:
            raise ValueError("Passed era should be either 'md', 'rd' or 'ded' ")

    @staticmethod
    @u.quantity_input(t=u.s)
    def scalar_factor_derivative(t, era="md", tuning_param=1.0):
        """
        Derivative of acceleration of the universe in cosmological models of Robertson Walker
        Flat Universe.

        Parameters
        ----------
        era : string
            Can be chosen from 'md' (Matter Dominant),
            'rd' (Radiation Dominant) and 'ded' (Dark Energy Dominant)
        t : ~astropy.units.s
            Time for the event
        tuning_param : float, optional
            Unit scaling factor, defaults to 1

        Returns
        -------
        float
            Value of derivative of scalar factor at time t.

        Raises
        ------
        ValueError : If era is not 'md' , 'rd', and 'ded'.

        """
        T = t.to(u.s).value
        if era == "md":
            return (2 / 3) * tuning_param * (T ** (-1 / 3))
        elif era == "rd":
            return 0.5 * tuning_param * (T ** (-0.5))
        elif era == "ded":
            hubble_const = (constant.Cosmo_Const / 3) ** 0.5
            val = hubble_const.value * (np.e ** (hubble_const.value * T))
            return tuning_param * val
        else:
            raise ValueError("Passed era should be either 'md', 'rd' or 'ded' ")
