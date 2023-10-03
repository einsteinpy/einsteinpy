import numpy as np
import sympy
from sympy import permutedims
from sympy.combinatorics import Permutation
from sympy.functions.special.tensor_functions import LeviCivita

from .christoffel import ChristoffelSymbols
from .einstein import EinsteinTensor
from .helpers import _change_name
from .levicivita import LeviCivitaAlternatingTensor
from .metric import MetricTensor
from .ricci import RicciScalar, RicciTensor
from .riemann import RiemannCurvatureTensor
from .stress_energy_momentum import StressEnergyMomentumTensor
from .tensor import BaseRelativityTensor, _change_config, tensor_product
from .vector import GenericVector
from .weyl import BelRobinsonTensor, WeylTensor

# from sympy import simplify, tensorcontraction, tensorproduct


class GenericSpacetime:
    """
    Class that encapsulates the metric and its associated tensors.
    All tensors needed during calculation are saved and accessible through properties.
    The properties can also be set, in case approximations are to be made along the way.
    """

    def __init__(
        self,
        metric,
        chris=None,
        riemann=None,
        ricci=None,
        einstein=None,
        sem_tensor=None,
        name="GenericSpacetime",
    ):
        """
        Constructor for the Generic Spacetime class. Only the MetricTensor is necessary,
        the other tensors can be provided if precomputed.

         Parameters
        ----------
        metric : ~einsteinpy.symbolic.metric.MetricTensor
            The metric tensor describing the spacetime

        chris : ~einsteinpy.symbolic.christoffel.ChristoffelSymbols  (optional)
            The associated Christoffel symbols
        riemann : ~einsteinpy.symbolic.riemann.RiemannCurvatureTensor  (optional)
            The associated Riemann curvature tensor
        ricci: : ~einsteinpy.symbolic.ricci.RicciTensor (optional)
            The associated Ricci tensor
        einstein : ~einsteinpy.symbolic.einstein.EinsteinTensor (optional)
            The associated Einstein tensor
        sem_tensor : ~einsteinpy.symbolic.stress_energy_momentum.StressEnergyMomentumTensor (optional)
            The associated SEM tensor
        name : String
            The name for the spacetime, default is "GenericSpacetime"
        """

        self._metric = metric
        self._chris = chris
        self._riemann = riemann
        self._ricci = ricci
        self._ricci_s = None
        self._einstein = einstein
        self._sem = sem_tensor
        self._levi_civita = None
        self._weyl = None
        self._belrob = None

    @property
    def EinsteinTensor(self):
        """
        Propery that returns the EinsteinTensor of the spacetime. Is calculated for the first time
        if not provided.

        Returns:
        -------
            _einstein : ~einsteinpy.symbolic.einstein.EinsteinTensor
                The Einstein tensor of the spacetime
        """
        if self._einstein is None:
            self._einstein = EinsteinTensor.from_ricci(
                self.RicciTensor, self.RicciScalar
            )
        return self._einstein

    @EinsteinTensor.setter
    def EinsteinTensor(self, value):
        self._einstein = value

    @property
    def RicciScalar(self):
        """
        Propery that returns the RicciScalar of the spacetime. Is calculated for the first time
        if not provided.

        Returns:
        -------
            _ricci_s : einsteinpy.symbolic.ricci.RicciScalar
                The Ricci scalar of the spacetime
        """
        if self._ricci_s is None:
            self._ricci_s = RicciScalar.from_riccitensor(self.RicciTensor)
        return self._ricci_s

    @RicciScalar.setter
    def RicciScalar(self, value):
        self._ricci_s = value

    @property
    def RicciTensor(self):
        """
        Propery that returns the RicciTensor of the spacetime. Is calculated for the first time
        if not provided.

        Returns:
        -------
            _ricci : ~einsteinpy.symbolic.ricci.RicciTensor
                The Ricci tensor of the spacetime
        """
        if self._ricci is None:
            self._ricci = RicciTensor.from_riemann(self.RiemannTensor)
        return self._ricci

    @RicciTensor.setter
    def RicciTensor(self, value):
        self._ricci = value

    @property
    def RiemannTensor(self):
        """
        Propery that returns the Riemann tensor of the spacetime. Is calculated for the first time
        if not provided.

        Returns:
        -------
            _riemann : ~einsteinpy.symbolic.riemann.RiemannTensor
                The Riemann tensor of the spacetime
        """
        if self._riemann is None:
            self._riemann = RiemannCurvatureTensor.from_christoffels(
                self.ChristoffelSymbols
            )
        return self._riemann

    @RiemannTensor.setter
    def RiemannTensor(self, value):
        self._riemann = value

    @property
    def WeylTensor(self):
        """
        Propery that returns the Weyl Tensor of the spacetime. Is calculated for the first time
        if not provided.

        Returns:
        -------
            _weyl : ~einsteinpy.symbolic.weyl.WeylTensor
                The Weyl tensor of the spacetime
        """
        if self._weyl is None:
            self._weyl = WeylTensor.from_tensors(
                self.Metric, self.RiemannTensor, self.RicciTensor, self.RicciScalar
            )
            self._weyl._parent_spacetime = self
        return self._weyl

    @WeylTensor.setter
    def WeylTensor(self, value):
        self._weyl = value

    @property
    def BelRobinsonTensor(self):
        """
        Propery that returns the BelRobinsonTensor of the spacetime.

        Returns:
        -------
            _belrob : ~einsteinpy.symbolic.belrobinson.BelRobinsonTensor
                The BelRobinsonTensor of the spacetime
        """
        if self._belrob is None:
            self._belrob = BelRobinsonTensor.from_weyl(self.WeylTensor)
        return self._belrob

    @BelRobinsonTensor.setter
    def BelRobinsonTensor(self, value):
        self._belrob = value

    @property
    def ChristoffelSymbols(self):
        """
        Propery that returns the Christoffel symbols of the spacetime. Is calculated for the first time
        if not provided.

        Returns:
        -------
            _chris : ~einsteinpy.symbolic.christoffel.ChristoffelSymbols
                The Christoffel symbols of the spacetime
        """
        if self._chris is None:
            self._chris = ChristoffelSymbols.from_metric(self.Metric)
        return self._chris

    @ChristoffelSymbols.setter
    def ChristoffelSymbols(self, value):
        self._chris = value

    @property
    def LeviCivitaTensor(self):
        """
        Propery that returns the alternating Levi-Civita tensor of the spacetime.

        Returns:
        -------
            _levi_civita : ~einsteinpy.symbolic.levicivita.LeviCivitaAlternatingTensor
                The alternating Levi-Civita tensor of the spacetime
        """
        if self._levi_civita is None:
            self._levi_civita = LeviCivitaAlternatingTensor.from_metric(self.Metric)
        return self._levi_civita

    @LeviCivitaTensor.setter
    def LeviCivitaTensor(self, value):
        self._levi_civita = value

    @property
    def SEMTensor(self):
        """
        Propery that returns the Stress-Energy-Momentum tensor of the spacetime. Is calculated for the first time
        if not provided.

        Returns:
        -------
            _sem : ~einsteinpy.symbolic.stress_energy_momentum.StressEnergyMomentumTensor
                The Stress-Energy-Momentum tensor of the spacetime
        """
        if self._sem is None:
            self._sem = StressEnergyMomentumTensor.from_einstein(self.EinsteinTensor)
        return self._sem

    @SEMTensor.setter
    def SEMTensor(self, value):
        self._sem = value

    @property
    def Metric(self):
        """
        Propery that returns the metric tensor of the spacetime.

        Returns:
        -------
            _einstein : ~einsteinpy.symbolic.metric.MetricTensor
                The metric tensor of the spacetime
        """
        return self._metric

    def geodesic_equation(self, wline, apar):
        """
        Returns the geodesic equation for a given wline x(apar) with the affine parameter apar

        Parameters
        --------
            wline : ~einsteinpy.symbolic.tensor.BaseRelativityTensor or GenericVector
                The 4d vector to describe the worldline
            apar : ~sympy.core.symbol.Symbol
                The symbol to describe the affine parameter to describe the worldline

        Returns
        ------
            ~sympy.core.relational.Equality
                The geodesic equation
        """
        chris = self.ChristoffelSymbols.change_config("ull")
        wline = wline.change_config("u")
        dx_ds = GenericVector(wline.tensor().diff(apar), syms=wline.syms, config="u")
        rhs = tensor_product(tensor_product(chris, dx_ds, 1, 0), dx_ds, 1, 0)
        d2x_ds2 = dx_ds.tensor().diff(apar)
        return sympy.Eq(d2x_ds2, -rhs.tensor())

    def covariant_derivative(self, T):
        """
        Calculates the covariant derivative for a given tensor or scalar T.

        Parameters
        -------
            T : ~einsteinpy.symbolic.tensor.BaseRelativityTensor or ~sympy.core.expr.Expr
                The tensor with arbitrary configuration, or a scalar.

        Returns
        -----
            Td : ~einsteinpy.symbolic.tensor.BaseRelativityTensor
                The covariant derivative of the given tensor or scalar

        Notes
        -----
            Tested for scalars, tensors of rank 1 and 2
        """
        chris = self.ChristoffelSymbols
        if not chris.config == "ull":
            chris = chris.change_config(newconfig="ull")

        syms = self.Metric.symbols()

        Tdel = []
        try:
            for s in syms:
                Tdel.append(T.tensor().diff(s))
            Td = BaseRelativityTensor(
                Tdel,
                syms=syms,
                config="l" + T.config,
                parent_metric=self.Metric,
                parent_spacetime=self,
            )

            for i in range(T.order):
                p = Permutation(size=T.order + 1)
                if T.config[i] == "u":
                    for j in range(i + 1):
                        p = Permutation(j, j + 1, size=T.order + 1) * p
                    Td.arr += permutedims(tensor_product(chris, T, 2, i).arr, p)
                if T.config[i] == "l":
                    p = Permutation(1, i + 1, size=T.order + 1) if i > 0 else p
                    Td.arr -= permutedims(tensor_product(chris, T, 0, i).arr, p)

            Td.simplify()
        except AttributeError:  # In case T is a scalar
            for s in syms:
                Tdel.append(T.diff(s))
            Td = BaseRelativityTensor(
                Tdel,
                syms=syms,
                config="l",
                parent_metric=self.Metric,
                parent_spacetime=self,
            )

        return Td
