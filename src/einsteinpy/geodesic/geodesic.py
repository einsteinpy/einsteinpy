import warnings

import numpy as np

from einsteinpy.integrators import GeodesicIntegrator

from .utils import _P, _kerr, _kerrnewman, _sch


class Geodesic:
    """
    Base Class for defining Geodesics
    Working in Geometrized Units (M-Units),
    with :math:`c = G = M = k_e = 1`

    """

    def __init__(
        self,
        metric,
        metric_params,
        position,
        momentum,
        time_like=True,
        return_cartesian=True,
        **kwargs,
    ):
        """
        Constructor

        Parameters
        ----------
        metric : str
            Name of the metric. Currently, these metrics are supported:
            1. Schwarzschild
            2. Kerr
            3. KerrNewman
        metric_params : array_like
            Tuple of parameters to pass to the metric
            E.g., ``(a,)`` for Kerr
        position : array_like
            3-Position
            4-Position is initialized by taking ``t = 0.0``
        momentum : array_like
            3-Momentum
            4-Momentum is calculated automatically,
            considering the value of ``time_like``
        time_like : bool, optional
            Determines type of Geodesic
            ``True`` for Time-like geodesics
            ``False`` for Null-like geodesics
            Defaults to ``True``
        return_cartesian : bool, optional
            Whether to return calculated positions in Cartesian Coordinates
            This only affects the coordinates. Momenta are dimensionless
            quantities, and are returned in Spherical Polar Coordinates.
            Defaults to ``True``
        kwargs : dict
            Keyword parameters for the Geodesic Integrator
            See 'Other Parameters' below.

        Other Parameters
        ----------------
        steps : int
            Number of integration steps
            Defaults to ``50``
        delta : float
            Initial integration step-size
            Defaults to ``0.5``
        rtol : float
            Relative Tolerance
            Defaults to ``1e-2``
        atol : float
            Absolute Tolerance
            Defaults to ``1e-2``
        order : int
            Integration Order
            Defaults to ``2``
        omega : float
            Coupling between Hamiltonian Flows
            Smaller values imply smaller integration error, but too
            small values can make the equation of motion non-integrable.
            For non-capture trajectories, ``omega = 1.0`` is recommended.
            For trajectories, that either lead to a capture or a grazing
            geodesic, a decreased value of ``0.01`` or less is recommended.
            Defaults to ``1.0``
        suppress_warnings : bool
            Whether to suppress warnings during simulation
            Warnings are shown for every step, where numerical errors
            exceed specified tolerance (controlled by ``rtol`` and ``atol``)
            Defaults to ``False``

        """
        # Contravariant Metrics, defined so far
        _METRICS = {
            "Schwarzschild": _sch,
            "Kerr": _kerr,
            "KerrNewman": _kerrnewman,
        }

        if metric not in _METRICS:
            if not callable(metric):
                raise NotImplementedError(
                    f"'{metric}' is unsupported. Currently, these metrics are supported:\
                    \n1. Schwarzschild\n2. Kerr\n3. KerrNewman\n4. Any user-supplied callable metric function"
                )

            self.metric_name = metric.__name__
            self.metric = metric
        else:
            self.metric_name = metric
            self.metric = _METRICS[metric]

        self.metric_params = metric_params
        if metric == "Schwarzschild":
            self.metric_params = (0.0,)
        self.position = np.array([0.0, *position])
        self.momentum = _P(
            self.metric, metric_params, self.position, momentum, time_like
        )
        self.time_like = time_like

        self.kind = "Time-like" if time_like else "Null-like"
        self.coords = "Cartesian" if return_cartesian else "Spherical Polar"

        self._trajectory = self.calculate_trajectory(**kwargs)

    def __repr__(self):
        return f"""Geodesic Object:(\n\
            Type : ({self.kind}),\n\
            Metric : ({self.metric_name}),\n\
            Metric Parameters : ({self.metric_params}),\n\
            Initial 4-Position : ({self.position}),\n\
            Initial 4-Momentum : ({self.momentum}),\n\
            Trajectory = (\n\
                {self.trajectory}\n\
            ),\n\
            Output Position Coordinate System = ({self.coords})\n\
        ))"""

    def __str__(self):
        return self.__repr__()

    @property
    def trajectory(self):
        """
        Returns the trajectory of the test particle

        """
        return self._trajectory

    def calculate_trajectory(self, **kwargs):
        """
        Calculate trajectory in spacetime

        Parameters
        ----------
        kwargs : dict
            Keyword parameters for the Geodesic Integrator
            See 'Other Parameters' below.

        Returns
        -------
        ~numpy.ndarray
            N-element numpy array, containing step count
        ~numpy.ndarray
            Shape-(N, 8) numpy array, containing
            (4-Position, 4-Momentum) for each step

        Other Parameters
        ----------------
        steps : int
            Number of integration steps
            Defaults to ``50``
        delta : float
            Initial integration step-size
            Defaults to ``0.5``
        rtol : float
            Relative Tolerance
            Defaults to ``1e-2``
        atol : float
            Absolute Tolerance
            Defaults to ``1e-2``
        order : int
            Integration Order
            Defaults to ``2``
        omega : float
            Coupling between Hamiltonian Flows
            Smaller values imply smaller integration error, but too
            small values can make the equation of motion non-integrable.
            For non-capture trajectories, ``omega = 1.0`` is recommended.
            For trajectories, that either lead to a capture or a grazing
            geodesic, a decreased value of ``0.01`` or less is recommended.
            Defaults to ``1.0``
        suppress_warnings : bool
            Whether to suppress warnings during simulation
            Warnings are shown for every step, where numerical errors
            exceed specified tolerance (controlled by ``rtol`` and ``atol``)
            Defaults to ``False``

        """
        g, g_prms = self.metric, self.metric_params
        q0, p0 = self.position, self.momentum
        tl = self.time_like
        N = kwargs.get("steps", 50)
        dl = kwargs.get("delta", 0.5)
        rtol = kwargs.get("rtol", 1e-2)
        atol = kwargs.get("atol", 1e-2)
        order = kwargs.get("order", 2)
        omega = kwargs.get("omega", 1.0)
        sw = kwargs.get("suppress_warnings", False)
        steps = np.arange(N)

        geodint = GeodesicIntegrator(
            metric=g,
            metric_params=g_prms,
            q0=q0,
            p0=p0,
            time_like=tl,
            steps=N,
            delta=dl,
            rtol=rtol,
            atol=atol,
            order=order,
            omega=omega,
            suppress_warnings=sw,
        )

        for i in steps:
            geodint.step()

        vecs = np.array(geodint.results, dtype=float)

        q1 = vecs[:, 0]
        p1 = vecs[:, 1]
        results = np.hstack((q1, p1))
        # Ignoring
        # q2 = vecs[:, 2]
        # p2 = vecs[:, 3]

        if self.coords == "Cartesian":
            # Converting to Cartesian from Spherical Polar Coordinates
            # Note that momenta cannot be converted this way,
            # due to ambiguities in the signs of v_r and v_th (velocities)
            t, r, th, ph = q1.T
            pt, pr, pth, pph = p1.T
            x = r * np.sin(th) * np.cos(ph)
            y = r * np.sin(th) * np.sin(ph)
            z = r * np.cos(th)

            cart_results = np.vstack((t, x, y, z, pt, pr, pth, pph)).T

            return steps, cart_results

        return steps, results


class Nulllike(Geodesic):
    """
    Class for defining Null-like Geodesics

    """

    def __init__(
        self, metric, metric_params, position, momentum, return_cartesian=True, **kwargs
    ):
        """
        Constructor

        Parameters
        ----------
        metric : str
            Name of the metric. Currently, these metrics are supported:
            1. Schwarzschild
            2. Kerr
            3. KerrNewman
        metric_params : array_like
            Tuple of parameters to pass to the metric
            E.g., ``(a,)`` for Kerr
        position : array_like
            3-Position
            4-Position is initialized by taking ``t = 0.0``
        momentum : array_like
            3-Momentum
            4-Momentum is calculated automatically,
            considering the value of ``time_like``
        return_cartesian : bool, optional
            Whether to return calculated positions in Cartesian Coordinates
            This only affects the coordinates. The momenta dimensionless
            quantities, and are returned in Spherical Polar Coordinates.
            Defaults to ``True``
        kwargs : dict
            Keyword parameters for the Geodesic Integrator
            See 'Other Parameters' below.

        Other Parameters
        ----------------
        steps : int
            Number of integration steps
            Defaults to ``50``
        delta : float
            Initial integration step-size
            Defaults to ``0.5``
        rtol : float
            Relative Tolerance
            Defaults to ``1e-2``
        atol : float
            Absolute Tolerance
            Defaults to ``1e-2``
        order : int
            Integration Order
            Defaults to ``2``
        omega : float
            Coupling between Hamiltonian Flows
            Smaller values imply smaller integration error, but too
            small values can make the equation of motion non-integrable.
            For non-capture trajectories, ``omega = 1.0`` is recommended.
            For trajectories, that either lead to a capture or a grazing
            geodesic, a decreased value of ``0.01`` or less is recommended.
            Defaults to ``1.0``
        suppress_warnings : bool
            Whether to suppress warnings during simulation
            Warnings are shown for every step, where numerical errors
            exceed specified tolerance (controlled by ``rtol`` and ``atol``)
            Defaults to ``False``

        """
        super().__init__(
            metric=metric,
            metric_params=metric_params,
            position=position,
            momentum=momentum,
            time_like=False,
            return_cartesian=return_cartesian,
            **kwargs,
        )


class Timelike(Geodesic):
    """
    Class for defining Time-like Geodesics

    """

    def __init__(
        self, metric, metric_params, position, momentum, return_cartesian=True, **kwargs
    ):
        """
        Constructor

        Parameters
        ----------
        metric : str
            Name of the metric. Currently, these metrics are supported:
            1. Schwarzschild
            2. Kerr
            3. KerrNewman
        metric_params : array_like
            Tuple of parameters to pass to the metric
            E.g., ``(a,)`` for Kerr
        position : array_like
            3-Position
            4-Position is initialized by taking ``t = 0.0``
        momentum : array_like
            3-Momentum
            4-Momentum is calculated automatically,
            considering the value of ``time_like``
        return_cartesian : bool, optional
            Whether to return calculated positions in Cartesian Coordinates
            This only affects the coordinates. The momenta dimensionless
            quantities, and are returned in Spherical Polar Coordinates.
            Defaults to ``True``
        kwargs : dict
            Keyword parameters for the Geodesic Integrator
            See 'Other Parameters' below.

        Other Parameters
        ----------------
        steps : int
            Number of integration steps
            Defaults to ``50``
        delta : float
            Initial integration step-size
            Defaults to ``0.5``
        rtol : float
            Relative Tolerance
            Defaults to ``1e-2``
        atol : float
            Absolute Tolerance
            Defaults to ``1e-2``
        order : int
            Integration Order
            Defaults to ``2``
        omega : float
            Coupling between Hamiltonian Flows
            Smaller values imply smaller integration error, but too
            small values can make the equation of motion non-integrable.
            For non-capture trajectories, ``omega = 1.0`` is recommended.
            For trajectories, that either lead to a capture or a grazing
            geodesic, a decreased value of ``0.01`` or less is recommended.
            Defaults to ``1.0``
        suppress_warnings : bool
            Whether to suppress warnings during simulation
            Warnings are shown for every step, where numerical errors
            exceed specified tolerance (controlled by ``rtol`` and ``atol``)
            Defaults to ``False``

        """
        super().__init__(
            metric=metric,
            metric_params=metric_params,
            position=position,
            momentum=momentum,
            time_like=True,
            return_cartesian=return_cartesian,
            **kwargs,
        )
