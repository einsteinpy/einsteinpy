import warnings

import numpy as np

from .utils import _Z, _flow_A, _flow_B, _flow_mixed


class GeodesicIntegrator:
    """
    Geodesic Integrator, based on [1]_.
    This module uses Forward Mode Automatic Differentiation
    to calculate metric derivatives to machine precision
    leading to stable simulations.

    References
    ----------
    .. [1] Christian, Pierre and Chan, Chi-Kwan;
        "FANTASY: User-Friendly Symplectic Geodesic Integrator
        for Arbitrary Metrics with Automatic Differentiation";
        `2021 ApJ 909 67 <https://doi.org/10.3847/1538-4357/abdc28>`__

    """

    # TODO: Update arXiv attributions to ApJ (See #572)
    def __init__(
        self,
        metric,
        metric_params,
        q0,
        p0,
        time_like=True,
        steps=100,
        delta=0.5,
        rtol=1e-2,
        atol=1e-2,
        order=2,
        omega=1.0,
        suppress_warnings=False,
    ):
        """
        Constructor

        Parameters
        ----------
        metric : callable
            Metric Function. Currently, these metrics are supported:
            1. Schwarzschild
            2. Kerr
            3. KerrNewman
        metric_params : array_like
            Tuple of parameters to pass to the metric
            E.g., ``(a,)`` for Kerr
        q0 : array_like
            Initial 4-Position
        p0 : array_like
            Initial 4-Momentum
        time_like : bool, optional
            Determines type of Geodesic
            ``True`` for Time-like geodesics
            ``False`` for Null-like geodesics
            Defaults to ``True``
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

        Raises
        ------
        NotImplementedError
            If ``order`` is not in [2, 4, 6, 8]

        """
        ORDERS = {
            2: self._ord_2,
            4: self._ord_4,
            6: self._ord_6,
            8: self._ord_8,
        }
        self.metric = metric
        self.metric_params = metric_params
        self.q0 = q0
        self.p0 = p0
        self.time_like = time_like
        self.steps = steps
        self.delta = delta
        self.omega = omega
        if order not in ORDERS:
            raise NotImplementedError(
                f"Order {order} integrator has not been implemented."
            )
        self.order = order
        self.integrator = ORDERS[order]
        self.rtol = rtol
        self.atol = atol
        self.suppress_warnings = suppress_warnings

        self.step_num = 0
        self.res_list = [q0, p0, q0, p0]
        self.results = list()

    def __str__(self):
        return f"""{self.__class__.__name__}(\n\
                metric : {self.metric}\n\
                metric_params : {self.metric_params}\n\
                q0 : {self.q0},\n\
                p0 : {self.p0},\n\
                time_like : {self.time_like},\n\
                steps : {self.steps},\n\
                delta : {self.delta},\n\
                omega : {self.omega},\n\
                order : {self.order},\n\
                rtol : {self.rtol},\n\
                atol : {self.atol}\n\
                suppress_warnings : {self.suppress_warnings}
            )"""

    def __repr__(self):
        return self.__str__()

    def _ord_2(self, q1, p1, q2, p2, delta):
        """
        Order 2 Integration Scheme

        References
        ----------
        .. [1] Christian, Pierre and Chan, Chi-Kwan;
            "FANTASY : User-Friendly Symplectic Geodesic Integrator
            for Arbitrary Metrics with Automatic Differentiation";
            `2021 ApJ 909 67 <https://doi.org/10.3847/1538-4357/abdc28>`__

        """
        dl, omg = delta, self.omega
        g = self.metric
        g_prms = self.metric_params

        HA1 = np.array(
            [
                q1,
                _flow_A(g, g_prms, q1, p1, q2, p2, 0.5 * dl)[1],
                _flow_A(g, g_prms, q1, p1, q2, p2, 0.5 * dl)[0],
                p2,
            ]
        )
        HB1 = np.array(
            [
                _flow_B(g, g_prms, HA1[0], HA1[1], HA1[2], HA1[3], 0.5 * dl)[0],
                HA1[1],
                HA1[2],
                _flow_B(g, g_prms, HA1[0], HA1[1], HA1[2], HA1[3], 0.5 * dl)[1],
            ]
        )
        HC = _flow_mixed(HB1[0], HB1[1], HB1[2], HB1[3], dl, omg)
        HB2 = np.array(
            [
                _flow_B(g, g_prms, HC[0], HC[1], HC[2], HC[3], 0.5 * dl)[0],
                HC[1],
                HC[2],
                _flow_B(g, g_prms, HC[0], HC[1], HC[2], HC[3], 0.5 * dl)[1],
            ]
        )
        HA2 = np.array(
            [
                HB2[0],
                _flow_A(g, g_prms, HB2[0], HB2[1], HB2[2], HB2[3], 0.5 * dl)[1],
                _flow_A(g, g_prms, HB2[0], HB2[1], HB2[2], HB2[3], 0.5 * dl)[0],
                HB2[3],
            ]
        )

        return HA2

    def _ord_4(self, q1, p1, q2, p2, delta):
        """
        Order 4 Integration Scheme

        References
        ----------
        .. [1] Yoshida, Haruo,
            "Construction of higher order symplectic integrators";
             Physics Letters A, vol. 150, no. 5-7, pp. 262-268, 1990.
            `DOI: <https://doi.org/10.1016/0375-9601(90)90092-3>`__

        """
        dl = delta

        Z0, Z1 = _Z(self.order)
        step1 = self._ord_2(q1, p1, q2, p2, dl * Z1)
        step2 = self._ord_2(step1[0], step1[1], step1[2], step1[3], dl * Z0)
        step3 = self._ord_2(step2[0], step2[1], step2[2], step2[3], dl * Z1)

        return step3

    def _ord_6(self, q1, p1, q2, p2, delta):
        """
        Order 6 Integration Scheme

        References
        ----------
        .. [1] Yoshida, Haruo,
            "Construction of higher order symplectic integrators";
             Physics Letters A, vol. 150, no. 5-7, pp. 262-268, 1990.
            `DOI: <https://doi.org/10.1016/0375-9601(90)90092-3>`__

        """
        dl = delta

        Z0, Z1 = _Z(self.order)
        step1 = self._ord_4(q1, p1, q2, p2, dl * Z1)
        step2 = self._ord_4(step1[0], step1[1], step1[2], step1[3], dl * Z0)
        step3 = self._ord_4(step2[0], step2[1], step2[2], step2[3], dl * Z1)

        return step3

    def _ord_8(self, q1, p1, q2, p2, delta):
        """
        Order 8 Integration Scheme

        References
        ----------
        .. [1] Yoshida, Haruo,
            "Construction of higher order symplectic integrators";
             Physics Letters A, vol. 150, no. 5-7, pp. 262-268, 1990.
            `DOI: <https://doi.org/10.1016/0375-9601(90)90092-3>`__

        """
        dl = delta

        Z0, Z1 = _Z(self.order)
        step1 = self._ord_6(q1, p1, q2, p2, dl * Z1)
        step2 = self._ord_6(step1[0], step1[1], step1[2], step1[3], dl * Z0)
        step3 = self._ord_6(step2[0], step2[1], step2[2], step2[3], dl * Z1)

        return step3

    def step(self):
        """
        Advances integration by one step

        """
        rl = self.res_list

        arr = self.integrator(rl[0], rl[1], rl[2], rl[3], self.delta)

        self.res_list = arr
        self.step_num += 1

        # Stability check
        if not self.suppress_warnings:
            g = self.metric
            g_prms = self.metric_params

            q1 = arr[0]
            p1 = arr[1]
            # Ignoring
            # q_2 = arr[2]
            # p_2 = arr[3]

            const = -int(self.time_like)
            # g.p.p ~ -1 or 0 (const)
            if not np.allclose(
                g(q1, *g_prms) @ p1 @ p1, const, rtol=self.rtol, atol=self.atol
            ):
                warnings.warn(
                    f"Numerical error has exceeded specified tolerance at step = {self.step_num}.",
                    RuntimeWarning,
                )

        self.results.append(self.res_list)
