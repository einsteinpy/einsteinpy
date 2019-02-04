import numpy as np
import warnings

class RK4():
    """
    Class for Defining Runge-Kutta 4th Order ODE solving method
    """
    def __init__(self, fun, t0, y0, t_bound, stepsize):
        """
        Initialization

        Parameters
        ----------
        fun : Function
            Should accept t, y as parameters, and return same type as y
        t0 : Float
            Initial t
        y0 : Numpy.array or Float
            Initial y
        t_bound : Float
            Boundary time - the integration won’t continue beyond it. It also determines the direction of the integration.
        stepsize : Float
            Size of each increment in t

        """
        self.t = t0
        self.y = y0
        self.t_bound = t_bound
        self.step_size = stepsize
        self._fun = fun
        self._postive_direction = (self.t_bound >= t0)
    
    def step(self):
        """
        Updates the value of self.t and self.y
        """
        if (self.t >= self.t_bound and self._postive_direction) or (self.t <= self.t_bound and not self._postive_direction):
            warnings.warn('Out of bounds set by t_bound. ', RuntimeWarning)
            return
        k0 = self._fun(self.t, self.y)
        k1 = self._fun(self.t + (self.step_size/2.0), self.y + (self.step_size/2.0)*k0)
        k2 = self._fun(self.t + (self.step_size/2.0), self.y + (self.step_size/2.0)*k1)
        k3 = self._fun(self.t + self.step_size, self.y + (self.step_size)*k2)
        self.y = self.y + ((self.step_size/6.0) * (k0 + 2*k1 + 2*k2 + k3))
        self.t = self.t + self.step_size
