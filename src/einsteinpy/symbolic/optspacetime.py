import numpy as np
import sympy
from sympy import simplify, tensorcontraction, tensorproduct

from .gem import GravitoElectricTensor, GravitoMagneticTensor
from .optdecomposition import OPTDecompositionTensor, OPTMetric
from .spacetime import GenericSpacetime, LeviCivitaAlternatingTensor
from .stress_energy_momentum import StressEnergyMomentumTensor
from .tensor import BaseRelativityTensor, Tensor, _change_name, tensor_product
from .vector import GenericVector


class ProjectedLeviCivitaAlternatingTensor(OPTDecompositionTensor):
    def __init__(
        self,
        arr,
        nvec,
        syms,
        config="lll",
        parent_metric=None,
        parent_spacetime=None,
        name="ProjectedAlternatingTensor",
    ):
        """
        Constructor for the Projected alternating Levi-Civita tensor

        Parameters
        ----------
        arr : ~sympy.tensor.array.dense_ndim_array.ImmutableDenseNDimArray or list
            Sympy Array or multi-dimensional list containing Sympy Expressions
        nvec : ~einsteinpy.symbolic.BaseRelativityTensor
            The normal unit time like vector used in the 1+3 decomposition
        syms : tuple or list
            Tuple of crucial symbols denoting time-axis, 1st, 2nd, and 3rd axis (t,x1,x2,x3)
        config : str
            Configuration of contravariant and covariant indices in tensor. 'u' for upper and 'l' for lower indices. Defaults to 'll'.
        parent_metric : ~einsteinpy.symbolic.metric.MetricTensor or None
            Corresponding Metric for the Ricci Tensor.
            Defaults to None.
        name : str
            Name of the Tensor

        Raises
        ------
        TypeError
            Raised when arr is not a list or sympy Array
        TypeError
            syms is not a list or tuple
        ValueError
            config has more or less than 2 indices

        """

        super(ProjectedLeviCivitaAlternatingTensor, self).__init__(
            arr,
            nvec=nvec,
            syms=syms,
            config=config,
            parent_metric=parent_metric,
            parent_spacetime=parent_spacetime,
            name=name,
        )

    @classmethod
    def from_metric(cls, optmetric):
        """
        Instantiates the LeviCivitaAlternatingTensor from a given metric

        Parameters
        ------
            optmetric : ~einsteinpy.symbolic.optdecomposition.OPTMetric
                The metric of the 1+3 decomposition

        Returns
        -------
            ~einsteinpy.symbolic.optspacetime.ProjectedLeviCivitaAlternatingTensor
        """
        lcat = LeviCivitaAlternatingTensor.from_metric(optmetric)
        eps = tensor_product(lcat, optmetric.NormalVector.change_config("u"), 3, 0)
        return ProjectedLeviCivitaAlternatingTensor(
            sympy.Array(eps.tensor()),
            nvec=optmetric.NormalVector,
            syms=optmetric.syms,
            config="lll",
            parent_metric=optmetric,
        )

    @classmethod
    def from_spacetime(cls, optst):
        """
        Instantiates the LeviCivitaAlternatingTensor from a given 1+3 decomposition spacetime

        Parameters
        ------
            optmetric : ~einsteinpy.symbolic.optspacetime.OPTSpacetime
                The spacetime object of the 1+3 decomposition

        Returns
        -------
            ~einsteinpy.symbolic.optspacetime.ProjectedLeviCivitaAlternatingTensor
        """
        optmetric = optst.Metric
        lcat = optst.LeviCivitaTensor
        eps = tensor_product(lcat, optmetric.NormalVector.change_config("u"), 3, 0)
        return ProjectedLeviCivitaAlternatingTensor(
            sympy.Array(eps.tensor()),
            nvec=optmetric.NormalVector,
            syms=optmetric.syms,
            config="lll",
            parent_metric=optmetric,
            parent_spacetime=optst,
        )


class OPTSEMTensor(OPTDecompositionTensor):
    """
    Class describing the Stress-Energy-Momentum Tensor in a 1+3 decomposition spacetime.
    The unit timelike vector nvec gives separates the spacetime into a

    Attributes:
    _rho : The energy density relativ to nvec
    _p   : The pressure relative to nvec
    _q   : The energy flux relative to nvec
    _pi  : The anisotropic stress relative to nvec
    """

    def __init__(
        self,
        arr,
        nvec,
        syms,
        config="ll",
        parent_metric=None,
        parent_spacetime=None,
        variables=list(),
        functions=list(),
        name="GenericOPTStressEnergyMomentumTensor",
    ):
        """
        The constructor for the OPTSEMTensor class

        Parameters
        ----------
        arr : ~sympy.tensor.array.dense_ndim_array.ImmutableDenseNDimArray or list
            Sympy Array or multi-dimensional list containing Sympy Expressions
        nvec : ~einsteinpy.symbolic.BaseRelativityTensor
            The normal unit time like vector used in the 1+3 decomposition
        syms : tuple or list
            Tuple of crucial symbols denoting time-axis, 1st, 2nd, and 3rd axis (t,x1,x2,x3)
        config : str
            Configuration of contravariant and covariant indices in tensor. 'u' for upper and 'l' for lower indices. Defaults to 'll'.
        parent_metric : ~einsteinpy.symbolic.metric.MetricTensor or None
            Metric Tensor associated with this Tensor.
        parent_spacetime : ~einsteinpy.symbolic.spacetime.GenericSpacetime or None
            Spacetime object associated with this Tensor.
        variables : tuple or list or set
            List of symbols used in expressing the tensor,
            other than symbols associated with denoting the space-time axis.
            Calculates in real-time if left blank.
        functions : tuple or list or set
            List of symbolic functions used in epressing the tensor.
            Calculates in real-time if left blank.
        name : str
            Name of the Tensor


        """
        super(OPTSEMTensor, self).__init__(
            arr,
            nvec=nvec,
            syms=syms,
            config=config,
            parent_metric=parent_metric,
            parent_spacetime=parent_spacetime,
            variables=variables,
            functions=functions,
            name=name,
        )

        self._rho = None
        self._p = None
        self._q = None
        self._pi = None

    @classmethod
    def from_einstein(cls, einstein, nvec=None):
        """
        Construction method to construct the SEM tensor simply with the Einstein tensor

        Parameters
        ---------
            einstein : ~einsteinpy.symbolic.einstein.EinsteinTensor
                The Einstein tensor
            nvec : ~einsteinpy.symbolic.tensor.BaseRelativityTensor

        Returns
        -------
            ~einsteinpy.symbolic.optspacetime.OPTSEMTensor
                The associated SEM tensor
        """
        sem = StressEnergyMomentumTensor.from_einstein(einstein)
        return cls(
            sem.arr,
            nvec=einstein.parent_metric.NormalVector if nvec is None else nvec,
            syms=sem.syms,
            config=sem.config,
            parent_metric=sem.parent_metric,
        )

    @property
    def rho(self):
        """
        Returns the energy density relative to the 1+3 decomposition.
        For the SEM tensor T_ab and the unit timelike vector u in the decomposition, this is given by
            rho = T_ab u^a u^b

        Returns
        ------
            _rho : sympy.core.expr.Expr
                The relative energy density
        """
        if self._rho is None:
            u = self.NormalVector.change_config("u")
            self._rho = tensorcontraction(
                tensor_product(
                    tensor_product(self.change_config("ll"), u, 0, 0), u
                ).arr,
                (0, 1),
            )
        return self._rho

    @rho.setter
    def rho(self, value):
        self._rho = value

    @property
    def p(self):
        """
        Returns the pressure relative to the 1+3 decomposition.
        For the SEM tensor T_ab and the projector tensor h_ab, this is given by
            rho = 1/3 T_ab h^ab

        Returns
        ------
            _p : sympy.core.expr.Expr
                The relative pressure
        """
        if self._p is None:
            h = self.parent_metric.ProjectorTensor.change_config("uu")
            self._p = (
                tensorcontraction(
                    tensor_product(self.change_config("ll"), h, 0, 0).arr, (0, 1)
                )
                / 3
            )
        return self._p

    @p.setter
    def p(self, value):
        self._p = value

    @property
    def q(self):
        """
        Returns the energy flux relative to the 1+3 decomposition.
        For the SEM tensor T_ab,the unit timelike vector u and the projector tensor h_ab in the decomposition, this is given by
            q_a = -h_a^c T_cb u^b

        Returns
        ------
            _q : ~einsteinpy.symbolic.optdecomposition.OPTDecompositionTensor
                The relative energy flux
        """
        if self._q is None:
            u = self.NormalVector.change_config("u")
            h = self.parent_metric.ProjectorTensor.change_config("lu")
            self._q = OPTDecompositionTensor(
                -tensor_product(
                    tensor_product(h, self.change_config("ll"), 1, 0), u, 1, 0
                ).tensor(),
                nvec=self.NormalVector,
                syms=self.syms,
                config="l",
                parent_metric=self.parent_metric,
            )
        return self._q

    @q.setter
    def q(self, value):
        self._q = value

    @property
    def pi(self):
        """
        Returns the anisotropic stress relative to the 1+3 decomposition.
        For the SEM tensor T_ab this is given by
            pi_ab = T_<ab>
        where <> signifies the projected symmetric trace free projection

        Returns
        ------
            _pi : ~einsteinpy.symbolic.optdecomposition.OPTDecompositionTensor
                The relative anisotropic stress
        """
        if self._pi is None:
            st = self._parent_spacetime or OPTSpacetime(self.parent_metric)
            self._pi = st.pstf(self)
        return self._pi

    @pi.setter
    def pi(self, value):
        self._pi = value


class OPTSpacetime(GenericSpacetime):
    """
    A class that encapsulates the spacetime of the 1+3 decomposition

    """

    def __init__(
        self,
        optmetric,
        chris=None,
        riemann=None,
        ricci=None,
        einstein=None,
        sem_tensor=None,
        name="OPTSpacetime",
    ):
        """
        Constructor for the OPTSpacetime class. Only the OPTMetricTensor is necessary,
        the other tensors can be provided if precomputed.

         Parameters
        ----------
        metric : ~einsteinpy.symbolic.optdecomposition.OPTMetricTensor
            The metric tensor describing the spacetime and 1+3 decomposition

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
            The name for the spacetime, default is "OPTSpacetime"
        """

        super(OPTSpacetime, self).__init__(
            metric=optmetric,
            chris=chris,
            riemann=riemann,
            ricci=ricci,
            einstein=einstein,
            sem_tensor=sem_tensor,
            name=name,
        )

        self._proj_alt_tensor = None
        self._proj_alt_tensor_dot = None
        self._expansion_scalar = None
        self._shear_tensor = None
        self._vorticity_tensor = None
        self._vorticity_vector = None
        self._E_tensor = None
        self._H_tensor = None
        self._u_dot = None
        self._Du = None

    def ProjectedNormalDerivatives(self):
        """
        Returns the projected derivatives of the unit timelike vector u of the 1+3 decomposition
            udot_a = u^b  \\Nabla_b u_a
            D_a u_b = h_a^c h_b^d \\Nabla_c u_d

        Returns
        ---------
            u_dot : einsteinpy.symbolic.OPTDecompositionTensor
        """
        if self._u_dot is None or self._Du is None:
            self._u_dot, self._Du = self.projected_covariant_derivative(
                self.NormalVector.change_config("l")
            )
        return self._u_dot, self._Du

    @property
    def ProjectedAlternatingTensor(self):
        """
        Returns the projected alternating Levi-Civita Tensor relative to the 1+3 decomposition.

        Returns
        ------
            ~einsteinpy.symbolic.optdecomposition.ProjectedLeviCivitaAlternatingTensor
        """
        if self._proj_alt_tensor is None:
            self._proj_alt_tensor = ProjectedLeviCivitaAlternatingTensor.from_spacetime(
                self
            )
        return self._proj_alt_tensor

    @ProjectedAlternatingTensor.setter
    def ProjectedAlternatingTensor(self, value):
        self._proj_alt_tensor = value

    @property
    def ProjectedAlternatingTensorDot(self):
        """
        Returns the temporal derivative of the projected alternating Levi-Civita Tensor relative to the 1+3 decomposition.

        Returns
        ------
            ~einsteinpy.symbolic.optdecomposition.OPTDecompositionTensor
        """
        if self._proj_alt_tensor_dot is None:
            self._proj_alt_tensor_dot, _ = self.projected_covariant_derivative(
                self.ProjectedAlternatingTensor
            )
        return self._proj_alt_tensor_dot

    @ProjectedAlternatingTensorDot.setter
    def ProjectedAlternatingTensorDot(self, value):
        self._proj_alt_tensor_dot = value

    @property
    def NormalVector(self):
        """
        Returns the normal unit timelike vector of the 1+3 decomposition.

        Returns
        ------
            ~einsteinpy.symbolic.vector.GenericVector
        """
        return self.Metric.NormalVector

    @property
    def ProjectorTensor(self):
        """
        Returns the projector tensor of the 1+3 decomposition.

        Returns
        ------
            ~einsteinpy.symbolic.optdecomposition.OPTDecompositionTensor
        """
        return self.Metric.ProjectorTensor

    @property
    def SEMTensor(self):
        """
        Returns the Stress-Energy-Momentum Tensor of the 1+3 decomposition.

        Returns
        ------
            ~einsteinpy.symbolic.optspacetime.OPTSEMTensor
        """
        if self._sem is None:
            self._sem = OPTSEMTensor.from_einstein(
                self.EinsteinTensor, nvec=self.Metric.NormalVector
            )
        return self._sem

    @SEMTensor.setter
    def SEMTensor(self, value):
        self._sem = value

    @property
    def ExpansionScalar(self):
        """
        Returns the expansion scalar of the 1+3 decomposition.

        Returns
        ------
            ~sympy.core.expr.Expr
        """
        if self._expansion_scalar is None:
            u_dot, Du = self.ProjectedNormalDerivatives()
            self._expansion_scalar = (
                tensorcontraction(Du.change_config("ul").tensor(), (0, 1))
                if self._expansion_scalar is None
                else self._expansion_scalar
            )
        return self._expansion_scalar

    @ExpansionScalar.setter
    def ExpansionScalar(self, value):
        self._expansion_scalar = value

    @property
    def ShearTensor(self):
        """
        Returns the shear tensor of the 1+3 decomposition.

        Returns
        ------
            ~einsteinpy.symbolic.optdecomposition.OPTDecompositionTensor
        """
        if self._shear_tensor is None:
            u_dot, Du = self.ProjectedNormalDerivatives()
            self._shear_tensor = self.pstf(Du)
        return self._shear_tensor

    @ShearTensor.setter
    def ShearTensor(self, value):
        self._shear_tensor = value

    @property
    def VorticityTensor(self):
        """
        Returns the vorticity tensor of the 1+3 decomposition.

        Returns
        ------
            ~einsteinpy.symbolic.optdecomposition.OPTDecompositionTensor
        """
        if self._vorticity_tensor is None:
            u_dot, Du = self.ProjectedNormalDerivatives()
            self._vorticity_tensor = OPTDecompositionTensor(
                Du.antisymmetric_part().tensor(),
                nvec=self.NormalVector,
                syms=self.Metric.syms,
                config="ll",
                parent_metric=self.Metric,
            )
        return self._vorticity_tensor

    @VorticityTensor.setter
    def VorticityTensor(self, value):
        self._vorticity_tensor = value

    @property
    def VorticityVector(self):
        """
        Returns the vorticity vector of the 1+3 decomposition.

        Returns
        ------
            ~einsteinpy.symbolic.optdecomposition.OPTDecompositionTensor
        """
        if self._vorticity_vector is None:
            self._vorticity_vector = OPTDecompositionTensor(
                -tensorcontraction(
                    tensor_product(
                        self.ProjectedAlternatingTensor,
                        self.VorticityTensor.change_config("uu"),
                        1,
                        0,
                    ).arr,
                    (1, 2),
                )
                / 2,
                nvec=self.NormalVector,
                syms=self.Metric.syms,
                config="l",
                parent_metric=self.Metric,
            )

        return self._vorticity_vector

    @VorticityVector.setter
    def VorticityVector(self, value):
        self._vorticity_vector = value

    @property
    def GravitoElectricTensor(self):
        """
        Returns the gravitoelectric tensor of the 1+3 decomposition.

        Returns
        ------
            ~einsteinpy.symbolic.gem.GravitoElectricTensor
        """
        if self._E_tensor is None:
            self._E_tensor = GravitoElectricTensor.from_weyl(
                self.WeylTensor, self.NormalVector
            )
        return self._E_tensor

    @GravitoElectricTensor.setter
    def GravitoElectricTensor(self, value):
        self._E_tensor = value

    @property
    def GravitoMagneticTensor(self):
        """
        Returns the gravitomagnetic tensor of the 1+3 decomposition.

        Returns
        ------
            ~einsteinpy.symbolic.gem.GravitoMagneticTensor
        """
        if self._H_tensor is None:
            self._H_tensor = GravitoMagneticTensor.from_opt_spacetime(self)
        return self._H_tensor

    @GravitoMagneticTensor.setter
    def GravitoMagneticTensor(self, value):
        self._H_tensor = value

    def spatial_trace(self, S):
        """
        Calculates the spatial trace of the given rank 2 tensor
            S = S_ab h^ab

        Parameters
        ---------
            S : einsteinpy.symbolic.tensor.BaseRelativityTensor
                Tensor of rank 2

        Returns
        -------
            sympy.core.expr.Expr
                The spatial trace
        """
        h = self.ProjectorTensor.change_config("uu")
        r = tensor_product(h, S.change_config("ll"), 0, 0)
        r = tensorcontraction(r.arr, (0, 1))
        return r

    def projection(self, T):
        """
        Calculates the spatial projection of any tensor
            P T_ab.. = h_a^c h_b^d T_cd..

        Parameters
        ---------
            S : einsteinpy.symbolic.tensor.BaseRelativityTensor
                Arbitrary Tensor

        Returns
        -------
            einsteinpy.symbolic.optdecomposition.OPTDecompositionTensor
                The spatial projection of the tensor
        """
        h = self.ProjectorTensor
        for i in range(T.order):
            if T.config[i] == "l":
                h = h.change_config("lu")
            elif T.config[i] == "u":
                h = h.change_config("ll")
            T = tensor_product(h, T, 1, i)
        return OPTDecompositionTensor(
            T.arr,
            nvec=self.Metric.NormalVector,
            syms=T.syms,
            config=T.config,
            parent_metric=self.Metric,
            parent_spacetime=self,
        )

    def projected_vector_dual(self, S):
        """
        Calculates the projected vector dual of the given rank 2 tensor
            S_a = 1/2 eps_abc S^[bc]

        Parameters
        ---------
            S : einsteinpy.symbolic.tensor.BaseRelativityTensor
                Tensor of rank 2

        Returns
        -------
            einsteinpy.symbolic.optdecomposition.OPTDecompositionTensor
                The projected vector dual
        """
        S = S.change_config("uu")
        S = S.antisymmetric_part()
        eps = self.ProjectedAlternatingTensor
        r = tensor_product(eps, S, 1, 0)
        r = tensorcontraction(r.arr, (1, 2)) / 2
        return OPTDecompositionTensor(
            r,
            nvec=self.Metric.NormalVector,
            syms=S.syms,
            config="l",
            parent_metric=self.Metric,
            parent_spacetime=self,
        )

    def pstf(self, S):
        """
        Calculates the projected symmetric tracefree tensor of the given rank 2 tensor

        Parameters
        ---------
            S : einsteinpy.symbolic.tensor.BaseRelativityTensor
                Tensor of rank 2

        Returns
        -------
            einsteinpy.symbolic.optdecomposition.OPTDecompositionTensor
                The associated projected symmetric tracefree tensor
        """
        h = self.ProjectorTensor
        S = S.change_config("ll")
        h_lu = h.change_config("lu")
        l = tensor_product(h_lu, S, 1, 0)
        l = tensor_product(h_lu, l, 1, 1)
        l = l.symmetric_part().arr
        r = tensor_product(h.change_config("uu"), S, 0, 0)
        r = tensorcontraction(r.arr, (0, 1))
        r = r * h.change_config("ll").arr
        return OPTDecompositionTensor(
            l - r / 3,
            nvec=self.Metric.NormalVector,
            syms=S.syms,
            config="ll",
            parent_metric=self.Metric,
            parent_spacetime=self,
        )

    def vector_product(self, V, W):
        """
        Calculates the vector product that is defined on the 1+3 decomposition
        Defined for [rank 1, rank 1], [rank 1, rank 2], or [rank 2, rank 2]

        Parameters
        --------
            V : einsteinpy.symbolic.tensor.BaseRelativityTensor
                Tensor of rank 1 or 2
            W : einsteinpy.symbolic.tensor.BaseRelativityTensor
                Tensor of rank 1 or 2

        Returns
        ------
            out : einsteinpy.symbolic.optdecomposition.OPTDecompositionTensor
                The vector product
        """
        eps = self.ProjectedAlternatingTensor
        if V.order == 1:
            if W.order == 1:
                r = tensor_product(eps, V.change_config("u"), 1, 0)
                r = tensor_product(r, W.change_config("u"), 1, 0)
                r = r.arr
            elif W.order == 2:
                S = W.change_config("lu")
                r = tensor_product(eps, S, 0, 1)
                r = tensor_product(r, V.change_config("u"), 1, 0)
                r = r.symmetric_part()
                r = r.arr
        elif V.order == 2:
            if W.order == 1:
                return self.vector_product(W, V)  # TODO: Check -
            S = V.change_config("ul")
            Q = W.change_config("uu")
            r = tensor_product(eps, S, 1, 0)
            r = tensor_product(r, Q, 2, 0)
            r = tensorcontraction(r.arr, (1, 2))
        else:
            raise Exception("Higher orders not implemented")
        return OPTDecompositionTensor(
            r,
            nvec=V.NormalVector,
            syms=V.syms,
            config="ll" if V.order == 1 and W.order == 2 else "l",
            parent_metric=self.Metric,
            parent_spacetime=self,
        )

    def projected_covariant_derivative(self, T):
        """
        Calculates the projected covariant derivatives of a given tensor or scalar T
            Tdot_a.. = u^b \\Nabla_b T_a..
            D_b T_a.. = h_b^d h_a^c \\Nabla_d T_c..

        Parameters
        ---------
            T : einsteinpy.symbolic.tensor.BaseRelativityTensor or sympy.core.expr.Expr
                The scalar or tensor

        Returns:
            Tdot : einsteinpy.symbolic.optdecomposition.OPTDecompositionTensor
                The time derivative
            DT : einsteinpy.symbolic.optdecomposition.OPTDecompositionTensor
                The spatial derivative
        """
        NablaT = self.covariant_derivative(T)
        u = self.Metric.NormalVector.change_config("u")
        h_lu = self.ProjectorTensor.change_config("lu")
        h_ul = self.ProjectorTensor.change_config("ul")
        Tdot = tensor_product(u, NablaT, 0, 0)
        DT = tensor_product(h_lu, NablaT, 1, 0)
        for i in range(1, DT.order):
            if DT.config[i] == "l":
                DT = tensor_product(DT, h_lu, i, 1)
            else:
                DT = tensor_product(DT, h_ul, i, 1)

        return (
            OPTDecompositionTensor(
                Tdot.arr,
                nvec=self.NormalVector,
                syms=T.syms,
                config=Tdot.config,
                parent_metric=self.Metric,
                parent_spacetime=self,
            ),
            OPTDecompositionTensor(
                DT.arr,
                nvec=self.NormalVector,
                syms=T.syms,
                config=DT.config,
                parent_metric=self.Metric,
                parent_spacetime=self,
            ),
        )

    def div(self, S):
        """
        Calculates the spatial divergence of a given tensor S of rank 1 or 2

        Parameters
        --------
            S : einsteinpy.symbolic.tensor.BaseRelativityTensor
                Tensor of rank 1 or 2

        Returns
        ------
            divS : einsteinpy.symbolic.optdecomposition.OPTDecompositionTensor or sympy.core.expr.Expr
                The spatial divergence
        """
        _, DS = self.projected_covariant_derivative(S)
        config = "u" + DS.config[1:]
        DS = DS.change_config(config)
        if S.order == 1:
            return tensorcontraction(DS.tensor(), (0, 1))
        else:
            return OPTDecompositionTensor(
                tensorcontraction(DS.tensor(), (0, S.order)),
                nvec=self.NormalVector,
                syms=S.syms,
                config=config[1:-1],
                parent_metric=self.Metric,
                parent_spacetime=self,
            )

    def curl(self, S):
        """
        Calculates the spatial curl of a given tensor S of rank 1 or 2

        Parameters
        --------
            S : einsteinpy.symbolic.tensor.BaseRelativityTensor
                Tensor of rank 1 or 2

        Returns
        ------
            curlS : einsteinpy.symbolic.optdecomposition.OPTDecompositionTensor
                The spatial curl
        """
        if S.order == 1:
            S = S.change_config("u")
        elif S.order == 2:
            S = S.change_config("lu")
        else:
            raise ValueError("Higher Orders not implemented")
        _, DS = self.projected_covariant_derivative(S)
        config = "u" + DS.config[1:]
        DS = DS.change_config(config)
        if S.order == 1:
            return OPTDecompositionTensor(
                tensorcontraction(
                    tensor_product(self.ProjectedAlternatingTensor, DS, 1, 0).arr,
                    (1, 2),
                ),
                nvec=S.NormalVector,
                syms=S.syms,
                config="l",
                parent_metric=self.Metric,
                parent_spacetime=self,
            )
        elif S.order == 2:
            return OPTDecompositionTensor(
                tensorcontraction(
                    tensor_product(self.ProjectedAlternatingTensor, DS, 0, 0).arr,
                    (0, 3),
                ),
                nvec=S.NormalVector,
                syms=S.syms,
                config="ll",
                parent_metric=self.Metric,
                parent_spacetime=self,
            ).symmetric_part()

    def temporal_rotation(self, S):
        """
        Calculates the temporal rotation of a given tensor of rank 1 or 2

        Parameters
        --------
            S : einsteinpy.symbolic.tensor.BaseRelativityTensor
                Tensor of rank 1 or 2

        Returns
        ------
            udotS : einsteinpy.symbolic.optdecomposition.OPTDecompositionTensor
                The temporal rotation
        """
        u = self.Metric.NormalVector.change_config("u")
        eps_dot = self.ProjectedAlternatingTensorDot
        if S.order == 1:
            r = tensor_product(eps_dot, S.change_config("u"), 1, 0)
            r = tensor_product(u, r, 0, 1)
            return OPTDecompositionTensor(
                -r.arr,
                nvec=u,
                syms=S.syms,
                config="l",
                parent_metric=self.Metric,
                parent_spacetime=self,
            )
        elif S.order == 2:
            r = tensor_product(eps_dot, S.change_config("lu"), 1, 1)
            r = tensor_product(u, r, 0, 0).symmetric_part()
            return OPTDecompositionTensor(
                -r.arr,
                nvec=u,
                syms=S.syms,
                config="ll",
                parent_metric=self.Metric,
                parent_spacetime=self,
            )
        else:
            raise ValueError("Higher Orders not implemented")
