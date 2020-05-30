import warnings

import numpy as np
from scipy import integrate


class RK4naive:
    """
    Class for Defining Runge-Kutta 4th Order ODE solving method
    """

    def __init__(self, fun, t0, y0, t_bound, stepsize):
        """
        Initialization

        Parameters
        ----------
        fun : function
            Should accept t, y as parameters, and return same type as y
        t0 : float
            Initial t
        y0 : ~numpy.array or float
            Initial y
        t_bound : float
            Boundary time - the integration won’t continue beyond it. It also determines the direction of the integration.
        stepsize : float
            Size of each increment in t

        """
        self.t = t0
        self.y = y0
        self.t_bound = t_bound
        self.step_size = stepsize
        self._fun = fun
        self.direction = 1 * (self.t_bound >= t0) - 1 * (not self.t_bound < t0)

    def step(self):
        """
        Updates the value of self.t and self.y
        """
        if (self.t >= self.t_bound and self.direction == 1) or (
            self.t <= self.t_bound and self.direction == -1
        ):
            warnings.warn("Out of bounds set by t_bound. ", RuntimeWarning)
            return
        k0 = self._fun(self.t, self.y)
        k1 = self._fun(
            self.t + (self.step_size / 2.0), self.y + (self.step_size / 2.0) * k0
        )
        k2 = self._fun(
            self.t + (self.step_size / 2.0), self.y + (self.step_size / 2.0) * k1
        )
        k3 = self._fun(self.t + self.step_size, self.y + (self.step_size) * k2)
        self.y = self.y + ((self.step_size / 6.0) * (k0 + 2 * k1 + 2 * k2 + k3))
        self.t = self.t + self.step_size


class RK45(integrate.RK45):
    """
    This Class inherits ~scipy.integrate.RK45 Class
    """

    def __init__(self, fun, t0, y0, t_bound, stepsize, rtol=None, atol=None):
        """
        Initialization

        Parameters
        ----------
        fun : function
            Should accept t, y as parameters, and return same type as y
        t0 : float
            Initial t
        y0 : ~numpy.array or float
            Initial y
        t_bound : float
            Boundary time - the integration won’t continue beyond it. It also determines the direction of the integration.
        stepsize : float
            Size of each increment in t
        rtol : float
            Relative tolerance, defaults to 0.2*stepsize
        atol : float
            Absolute tolerance, defaults to rtol/0.8e3
        
        """
        vectorized = not isinstance(y0, float)
        self._t_bound = t_bound
        if rtol is None:
            rtol = 0.2 * stepsize
        if atol is None:
            atol = rtol / 0.8e3
        super(RK45, self).__init__(
            fun=fun,
            t0=t0,
            y0=y0,
            t_bound=self._t_bound,
            first_step=0.8 * stepsize,
            max_step=8.0 * stepsize,
            rtol=rtol,
            atol=atol,
            vectorized=vectorized,
        )

    def step(self):
        """
        Updates the value of self.t and self.y
        """

        try:
            super(RK45, self).step()
        except RuntimeError:
            warnings.warn(
                "Attempt to step on a failed or finished solver. (Invalid Value or out of bounds of t_bound)",
                RuntimeWarning,
            )


# DRAFT CHANGES/ADDITIONS -- STARTS HERE - ????
# Cannot make adaptive meshing an optional parameter
# in the RK45 class, above, as stepsize control is
# not available in scipy.integrate.
# In prep for use with KerrShadow in rays.shadow
class RK45Adaptive:  #
    """
    Implements RK45 integration scheme, with Adaptive Stepsize Control
    """

    def __init__(self, fun, t0, y0, t_bound, stepsize, rtol=None, atol=None):
        """
        Initialization

        Parameters
        ----------
        fun : function
            Should accept t, y as parameters, and return same type as y
        t0 : float
            Initial t (Affine parameter, Proper Time for timelike geodesics) - ????
        y0 : ~numpy.array or float
            Initial State Vector - ????
        t_bound : float
            Integration Bounds - the integration won’t continue beyond it. It also determines the direction of the integration.
        stepsize : float
            Size of initial step - ????
        
        # Default values have to chosen balancing accuracy and performance - ????
        rtol : float
            Relative tolerance, defaults to 0.2*stepsize
        atol : float
            Absolute tolerance, defaults to rtol/0.8e3
        
        """
        pass

    def _nextstepsize(self, x1, x2, x3, u1, u2, u3, prev_step):
        """
        Returns next stepsize, which depends on position (x), velocity (u) and previous stepsize
        Specific to BL Coordinates, for now - ????

        Based on RAPTOR (2018) - ????

        Parameters
        ----------
        x1 : 
            Value of X^1 Coordinate (r, in BL Coordinates)
        x2 : 
            Value of X^2 Coordinate (theta, in BL Coordinates)
        x3 : 
            Value of X^3 Coordinate (phi, in BL Coordinates)
        u1 : 
            Value of U^1 Coordinate (u_r, in BL Coordinates)
        u2 : 
            Value of U^2 Coordinate (u_theta, in BL Coordinates)
        u3 : 
            Value of U^3 Coordinate (u_phi, in BL Coordinates)
        prev_step : float
            Size of previous integration step
        
        Returns
        -------
        float :
            Stepsize for next integration step
            
        """
        # To protect against division by 0
        delta = 1e-20  # ???? Choose appropriate delta
        delta2 = delta ** 2

        dlx1 = prev_step / (np.abs(u1) + delta2)
        dlx2 = prev_step * np.min(x2, 1.0 - x2) / (np.abs(u2) + delta2)
        dlx3 = prev_step / (np.abs(u3) + delta2)

        idlx1 = 1.0 / (np.abs(dlx1) + delta2)
        idlx2 = 1.0 / (np.abs(dlx2) + delta2)
        idlx3 = 1.0 / (np.abs(dlx3) + delta2)

        next_step = -np.max([1.0 / (idlx1 + idlx2 + idlx3), 1e-12])
        return next_step

    def step(self):
        """
        Moves integration forward by one step
        """
<<<<<<< HEAD
        # Calls _nextstepsize() for getting the stepsize of next integration step
=======
        # Calls _nextstep() for getting the stepsize
>>>>>>> fc0ed0a... Draft PR - integrator | metric
        pass


# DRAFT CHANGES/ADDITIONS -- ENDS HERE - ????
