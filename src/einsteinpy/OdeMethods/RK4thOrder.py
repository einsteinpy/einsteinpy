import numpy as np

class RK4thOrder():
    """
    Class for Defining Runge-Kutta 4th Order ODE solving method
    """
    def __init__(self, fun, t0, y0, stepsize):
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
        stepsize : Float
            Size of each increment in t

        """
        self.t = t0
        self.y = y0
        self.h = stepsize
        self._fun = fun
    
    def step(self):
        """
        Updates the value of self.t and self.y
        """
        k0 = self._fun(self.t, self.y)
        k1 = self._fun(self.t + (self.h/2.0), self.y + (self.h/2.0)*k0)
        k2 = self._fun(self.t + (self.h/2.0), self.y + (self.h/2.0)*k1)
        k3 = self._fun(self.t + self.h, self.y + (self.h)*k2)
        self.y = self.y + ((self.h/6.0) * (k0 + 2*k1 + 2*k2 + k3))
        self.t = self.t + self.h