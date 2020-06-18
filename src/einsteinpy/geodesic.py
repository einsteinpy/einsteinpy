import warnings

import numpy as np

from einsteinpy.coordinates import BoyerLindquistConversion, SphericalConversion
from einsteinpy.integrators import RK45


class Geodesic:
    """
    Base Class for defining Geodesics

    """

    def __init__(self, metric, init_vec, end_lambda, step_size=1e-3):
        """
        Parameters
        ----------
        metric : ~einsteinpy.metric.*
            Metric, in which Geodesics are to be calculated
        init_vec : numpy.array
            Length-8 Vector containing Initial 4-Position and 4-Velocity, \
            in that order
        end_lambda : float
            Affine Parameter, Lambda, where iterations will stop
            Equivalent to Proper Time for Timelike Geodesics
        step_size : float, optional
            Size of each geodesic integration step
            Defaults to ``1e-3``

        """
        self.metric = metric
        self.init_vec = init_vec
        # Showing messages, mainly in cases, when calculation is lengthy
        print("Calculating geodesic...")
        self._trajectory = self.calculate_trajectory(
            end_lambda=end_lambda, OdeMethodKwargs={"stepsize": step_size}
        )[1]
        print("Done!")

    def __repr__(self):
        return f"Geodesic:\n\nMetric = ({self.metric}),\n\ninit_vec = ({self.init_vec}),\n\nTrajectory = ({self.trajectory})"

    def __str__(self):
        return f"Geodesic:\n\nMetric = ({self.metric}),\n\ninit_vec = ({self.init_vec}),\n\nTrajectory = ({self.trajectory})"

    @property
    def trajectory(self):
        return self._trajectory

    def calculate_trajectory(
        self,
        end_lambda=10.0,
        stop_on_singularity=True,
        OdeMethodKwargs={"stepsize": 1e-3},
        return_cartesian=True,
    ):
        """
        Calculate trajectory in spacetime, according to Geodesic Equations
        
        Parameters
        ----------
        end_lambda : float
            Affine Parameter, Lambda, where iterations will stop
            Equivalent to Proper Time for Timelike Geodesics
            Defaults to ``10``
        stop_on_singularity : bool
            Whether to stop further computation on reaching singularity
            Defaults to ``True``
        OdeMethodKwargs : dict
            Kwargs to be supplied to the ODESolver
            Dictionary with key 'stepsize' along with a float value is expected
            Defaults to ``{'stepsize': 1e-3}``
        return_cartesian : bool
            Whether to return calculated values in Cartesian Coordinates
            Defaults to ``True``

        Returns
        -------
        ~numpy.ndarray
            N-element array containing Lambda, where the geodesic equations were evaluated
        ~numpy.ndarray
            (n,8) shape array containing [x0, x1, x2, x3, v0, v1, v2, v3] for each Lambda

        """
        ODE = RK45(
            fun=self.metric.f_vec,
            t0=0.0,
            y0=self.init_vec,
            t_bound=end_lambda,
            **OdeMethodKwargs,
        )

        vecs = list()
        lambdas = list()
        crossed_event_horizon = False
        _scr = self.metric.sch_rad * 1.001

        while ODE.t < end_lambda:
            vecs.append(ODE.y)
            lambdas.append(ODE.t)
            ODE.step()
            if (not crossed_event_horizon) and (ODE.y[1] <= _scr):
                warnings.warn(
                    "Test particle has reached Schwarzchild Radius. ", RuntimeWarning
                )
                if stop_on_singularity:
                    break
                else:
                    crossed_event_horizon = True

        vecs, lambdas = np.array(vecs), np.array(lambdas)

        if not return_cartesian:
            return lambdas, vecs

        elif self.metric.coords == "S":
            cart_vecs = list()
            for v in vecs:
                si_vals = SphericalConversion(
                    v[1], v[2], v[3], v[5], v[6], v[7]
                ).convert_cartesian()
                cart_vecs.append(np.hstack((v[0], si_vals[:3], v[4], si_vals[3:])))

            return lambdas, np.array(cart_vecs)

        elif self.metric.coords == "BL":
            cart_vecs = list()
            for v in vecs:
                si_vals = BoyerLindquistConversion(
                    v[1], v[2], v[3], v[5], v[6], v[7]
                ).convert_cartesian()
                cart_vecs.append(np.hstack((v[0], si_vals[:3], v[4], si_vals[3:])))

            return lambdas, np.array(cart_vecs)

    def calculate_trajectory_iterator(
        self,
        stop_on_singularity=True,
        OdeMethodKwargs={"stepsize": 1e-3},
        return_cartesian=False,
    ):
        """
        Calculate trajectory in manifold according to geodesic equation
        Yields an iterator

        Parameters
        ----------
        stop_on_singularity : bool
            Whether to stop further computation on reaching singularity
            Defaults to ``True``
        OdeMethodKwargs : dict
            Kwargs to be supplied to the ODESolver, defaults to {'stepsize': 1e-3}
            Dictionary with key 'stepsize' along with an float value is expected.
        return_cartesian : bool
            Whether to return calculated values in Cartesian Coordinates
            Defaults to ``True``

        Yields
        ------
        float
            proper time
        ~numpy.ndarray
            array of [t, x1, x2, x3, velocity_of_time, v1, v2, v3] for each proper time(lambda).

        """
        ODE = RK45(
            fun=self.metric.f_vec,
            t0=0.0,
            y0=self.init_vec,
            t_bound=1e300,
            **OdeMethodKwargs,
        )

        crossed_event_horizon = False
        _scr = self.metric.sch_rad * 1.001

        while True:
            if not return_cartesian:
                yield ODE.t, ODE.y
            elif self.metric.coords == "S":
                v = ODE.y
                si_vals = SphericalConversion(
                    v[1], v[2], v[3], v[5], v[6], v[7]
                ).convert_cartesian()
                yield ODE.t, np.hstack((v[0], si_vals[:3], v[4], si_vals[3:]))
            elif self.metric.coords == "BL":
                v = ODE.y
                si_vals = BoyerLindquistConversion(
                    v[1], v[2], v[3], v[5], v[6], v[7]
                ).convert_cartesian()
                yield ODE.t, np.hstack((v[0], si_vals[:3], v[4], si_vals[3:]))
            ODE.step()
            if (not crossed_event_horizon) and (ODE.y[1] <= _scr):
                warnings.warn(
                    "Test particle has reached Schwarzchild Radius. ", RuntimeWarning
                )
                if stop_on_singularity:
                    break
                else:
                    crossed_event_horizon = True
