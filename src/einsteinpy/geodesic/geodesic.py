import warnings

import numpy as np

from einsteinpy.integrators import RK45


class Geodesic:
    """
    Base Class for defining Geodesics

    """

    # Add a `state` class variable, which will
    # store a function supplied by the user.
    # For usual cases, it will be handled
    # automatically by Timelike and NullGeodesic
    # classes. This variable is expected to be used
    # only for user-defined Geodesic classes. - ?????
    def __init__(
        self,
        time_like,
        metric,
        coords,
        end_lambda,
        step_size=1e-3,
        return_cartesian=True,
    ):
        """
        Parameters
        ----------
        time_like : bool
            Determines type of Geodesic
            ``True`` for Time-like geodesics
            ``False`` for Null-like geodesics
        metric : ~einsteinpy.metric.*
            Metric, in which Geodesics are to be calculated
        coords : ~einsteinpy.coordinates.differential.*
            Coordinate system, in which Metric is to be represented
        end_lambda : float
            Affine Parameter, Lambda, where iterations will stop
            Equivalent to Proper Time for Timelike Geodesics
        step_size : float, optional
            Size of each geodesic integration step
            Defaults to ``1e-3``
        return_cartesian : bool, optional
            Whether to return calculated values in Cartesian Coordinates
            Defaults to ``True``

        """
        self.time_like = time_like
        self.metric = metric
        self.coords = coords

        self._state = self._calculate_state()
        self._trajectory = self.calculate_trajectory(
            end_lambda=end_lambda,
            OdeMethodKwargs={"stepsize": step_size},
            return_cartesian=return_cartesian,
        )[1]

    def __repr__(self):
        kind = "Null-like"
        if self.time_like:
            kind = "Time-like"

        return f"Geodesic Object:\n\
            Type = ({kind})\n\
            Metric = ({self.metric}),\n\
            Initial State = ({self.state}),\n\
            Trajectory = ({self.trajectory})"

    def __str__(self):
        kind = "Null-like"
        if self.time_like:
            kind = "Time-like"

        return f"Geodesic Object:\n\
            Type = ({kind})\n\
            Metric = ({self.metric}),\n\
            Initial State = ({self.state}),\n\
            Trajectory = ({self.trajectory})"

    @property
    def state(self):
        """
        Returns the Initial State Vector of the Geodesic

        """
        return self._state

    @property
    def trajectory(self):
        """
        Returns the "Trajectory" of the Geodesic

        """
        return self._trajectory

    # Move to Timelike - ?????
    def _calculate_state(self):
        """
        Prepares and returns the Initial State Vector of the test particle

        Returns
        -------
        state : ~numpy.ndarray
            Initial State Vector of the test particle
            Length-8 Array

        Raises
        ------
        TypeError
            If there is a mismatch between the coordinates class of ``self.coords`` and \
            coordinate class, with which ``self.metric`` object has been instantiated

        """
        if self.coords.system != self.metric.coords.system:
            raise TypeError(
                "Coordinate System Mismatch between Metric object and supplied initial coordinates."
            )

        x4 = self.coords.position()
        v4 = self.coords.velocity(self.metric)

        state = np.hstack((x4, v4))

        return state

    def calculate_trajectory(
        self,
        end_lambda=10.0,
        OdeMethodKwargs={"stepsize": 1e-3},
        return_cartesian=True,
    ):
        """
        Calculate trajectory in spacetime, according to Geodesic Equations

        Parameters
        ----------
        end_lambda : float, optional
            Affine Parameter, Lambda, where iterations will stop
            Equivalent to Proper Time for Timelike Geodesics
            Defaults to ``10.0``
        OdeMethodKwargs : dict, optional
            Kwargs to be supplied to the ODESolver
            Dictionary with key 'stepsize' along with a float value is expected
            Defaults to ``{'stepsize': 1e-3}``
        return_cartesian : bool, optional
            Whether to return calculated values in Cartesian Coordinates
            Defaults to ``True``

        Returns
        -------
        ~numpy.ndarray
            N-element numpy array containing Lambda, where the geodesic equations were evaluated
        ~numpy.ndarray
            (n,8) shape numpy array containing [x0, x1, x2, x3, v0, v1, v2, v3] for each Lambda

        """
        ODE = RK45(
            fun=self.metric.f_vec,
            t0=0.0,
            y0=self.state,
            t_bound=end_lambda,
            **OdeMethodKwargs,
        )

        g = self.metric

        vecs = list()
        lambdas = list()
        _scr = g.sch_rad * 1.001

        while ODE.t < end_lambda:
            vecs.append(ODE.y)
            lambdas.append(ODE.t)
            ODE.step()

            if ODE.y[1] <= _scr:
                warnings.warn(
                    "Test particle has reached Schwarzchild Radius. ", RuntimeWarning
                )
                break

        vecs, lambdas = np.array(vecs), np.array(lambdas)

        if return_cartesian:
            # Getting the corresponding Conversion superclass
            # of the differential object
            conv = type(self.coords).__bases__[0]
            cart_vecs = list()
            for v in vecs:
                vals = conv(v[0], v[1], v[2], v[3], v[5], v[6], v[7]).convert_cartesian(
                    M=g.M, a=g.a
                )
                cart_vecs.append(np.hstack((vals[:4], v[4], vals[4:])))

            return lambdas, np.array(cart_vecs)

        return lambdas, vecs

    def calculate_trajectory_iterator(
        self, OdeMethodKwargs={"stepsize": 1e-3}, return_cartesian=True,
    ):
        """
        Calculate trajectory in manifold according to geodesic equation
        Yields an iterator

        Parameters
        ----------
        OdeMethodKwargs : dict, optional
            Kwargs to be supplied to the ODESolver
            Dictionary with key 'stepsize' along with a float value is expected
            Defaults to ``{'stepsize': 1e-3}``
        return_cartesian : bool, optional
            Whether to return calculated values in Cartesian Coordinates
            Defaults to ``True``

        Yields
        ------
        float
            Affine Parameter, Lambda, where the geodesic equations were evaluated
        ~numpy.ndarray
            Numpy array containing [x0, x1, x2, x3, v0, v1, v2, v3] for each Lambda

        """
        ODE = RK45(
            fun=self.metric.f_vec,
            t0=0.0,
            y0=self.state,
            t_bound=1e300,
            **OdeMethodKwargs,
        )

        g = self.metric

        _scr = g.sch_rad * 1.001

        while True:
            if return_cartesian:
                # Getting the corresponding Conversion superclass
                # of the differential object
                conv = type(self.coords).__bases__[0]
                v = ODE.y
                vals = conv(v[0], v[1], v[2], v[3], v[5], v[6], v[7]).convert_cartesian(
                    M=g.M, a=g.a
                )

                yield ODE.t, np.hstack((vals[:4], v[4], vals[4:]))

            else:
                yield ODE.t, ODE.y

            ODE.step()

            if ODE.y[1] <= _scr:
                warnings.warn(
                    "Test particle has reached Schwarzchild Radius. ", RuntimeWarning
                )
                break
