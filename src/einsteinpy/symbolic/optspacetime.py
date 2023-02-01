import numpy as np
import sympy
from sympy import simplify, tensorcontraction, tensorproduct

from .tensor import BaseRelativityTensor, Tensor, tensor_product, _change_name
from .vector import GenericVector
from .spacetime import GenericSpacetime, LeviCivitaAlternatingTensor
from .optdecomposition import OPTDecompositionTensor, OPTMetric
from .gem import GravitoElectricTensor, GravitoMagneticTensor
from .stress_energy_momentum import StressEnergyMomentumTensor


class ProjectedAlternatingTensor(OPTDecompositionTensor):

    def __init__(self, optmetric):
        lcat = LeviCivitaAlternatingTensor.from_metric(optmetric)
        eps = tensor_product(lcat, optmetric.normal_vector.change_config("u"), 3, 0)
        super(ProjectedAlternatingTensor, self).__init__(sympy.Array(eps.tensor()),  nvec=optmetric.normal_vector, syms=optmetric.syms, config="lll", parent_metric=optmetric)


    '''
    def GetDualTensor(self, T):

        T = T.change_config("llll")
        eps = self.change_config( "uull")

        dual = 1./2. * sympy.tensorproduct(T.arr, eps.arr)
        dual = sympy.simplify(sympy.tensorcontraction(sympy.tensorcontraction(dual, (2, 4)), (2, 4)))
        return BaseRelativityTensor(dual, syms=T.syms, config="llll", parent_metric=self.parent_metric)
    '''

class OPTSEMTensor(OPTDecompositionTensor):
    """

    """

    def __init__(
        self,
        arr,
        nvec,
        syms,
        config="ll",
        parent_metric=None,
        variables=list(),
        functions=list(),
        name="GenericOPTStressEnergyMomentumTensor",
    ):
        super(OPTSEMTensor, self).__init__(arr, nvec=nvec, syms=syms, config=config, parent_metric=parent_metric, variables=variables, functions=functions, name=name)

        self._rho = None
        self._p = None
        self._q = None
        self._pi = None

    @classmethod
    def from_einstein(cls, einstein, nvec=None):
        sem = StressEnergyMomentumTensor.from_einstein(einstein)
        return cls(sem.arr, nvec=einstein.parent_metric.normal_vector if nvec is None else nvec, syms=sem.syms, config=sem.config, parent_metric=sem.parent_metric)

    @property
    def rho(self):
        if self._rho is None:
            u = self.normal_vector.change_config("u")
            self._rho = tensorcontraction(tensor_product(tensor_product(self.change_config("ll"), u, 0, 0), u).arr, (0,1))
        return self._rho

    @rho.setter
    def rho(self, value):
        self._rho = value

    @property
    def p(self):
        if self._p is None:
            h = self.parent_metric.projector_tensor.change_config("uu")
            self._p = tensorcontraction(tensor_product(self.change_config("ll"), h, 0, 0).arr, (0,1))/3
        return self._p

    @p.setter
    def p(self, value):
        self._p = value

    @property
    def q(self):
        if self._q is None:
            u = self.normal_vector.change_config("u")
            h = self.parent_metric.projector_tensor.change_config("lu")
            self._q = OPTDecompositionTensor( -tensor_product(tensor_product(h, self.change_config("ll"), 1, 0), u, 1, 0).tensor(),
                                                nvec=self.normal_vector, syms=self.syms, config="l", parent_metric=self.parent_metric)
        return self._q

    @q.setter
    def q(self, value):
        self._q = value

    @property
    def pi(self):
        if self._pi is None:
            st = OPTSpacetime(self.parent_metric)
            self._pi = st.pstf(self)
        return self._pi

    @pi.setter
    def pi(self, value):
        self._pi = value

class OPTSpacetime(GenericSpacetime):

    def __init__( self,
                    optmetric,
                    chris=None,
                    riemann=None,
                    ricci=None,
                    einstein=None,
                    sem_tensor=None,
                    name = "OPTSpacetime"):

        super(OPTSpacetime, self).__init__(metric=optmetric, chris=chris, riemann=riemann, ricci=ricci, einstein=einstein, sem_tensor=sem_tensor, name=name)

        self._proj_alt_tensor = None
        self._proj_alt_tensor_dot = None
        self._acceleration_vector = None
        self._expansion_scalar = None
        self._shear_tensor = None
        self._vorticity_tensor = None
        self._vorticity_vector = None
        self._E_tensor = None
        self._H_tensor = None


    def _calculate_kinematic_quantities(self):
        u_down = self.NormalVector.change_config("l")
        u_dot, Du = self.projected_covariant_derivative(u_down)
        self._acceleration_vector = u_dot
        self._expansion_scalar = tensorcontraction(Du.change_config("ul").tensor(), (0,1))
        self._shear_tensor = self.pstf(Du)
        self._vorticity_tensor = Du.antisymmetric_part()
        self._vorticity_vector = OPTDecompositionTensor( -
                                    tensorcontraction(
                                            tensor_product(self.ProjectedAlternatingTensor, self._vorticity_tensor.change_config("uu"), 1, 0).arr,
                                                        (1,2))/2,
                                                        nvec=self.NormalVector, syms=self.Metric.syms, config="l", parent_metric=self.Metric)



    @property
    def ProjectedAlternatingTensor(self):
        if self._proj_alt_tensor is None:
            self._proj_alt_tensor = ProjectedAlternatingTensor(self.Metric)
        return self._proj_alt_tensor

    @ProjectedAlternatingTensor.setter
    def ProjectedAlternatingTensor(self, value):
        self._proj_alt_tensor = value

    @property
    def ProjectedAlternatingTensorDot(self):
        if self._proj_alt_tensor_dot is None:
            self._proj_alt_tensor_dot = tensor_product(self.Metric.normal_vector.change_config("u"), self.covariant_derivative(self.ProjectedAlternatingTensor), 0, 0)
        return self._proj_alt_tensor_dot

    @ProjectedAlternatingTensorDot.setter
    def ProjectedAlternatingTensorDot(self, value):
        self._proj_alt_tensor_dot = value

    @property
    def NormalVector(self):
        return self.Metric.normal_vector

    @property
    def ProjectorTensor(self):
        return self.Metric.projector_tensor

    @property
    def SEMTensor(self):
        if self._sem_tensor is None:
            self._sem_tensor = OPTSEMTensor.from_einstein(self.EinsteinTensor, nvec=self.Metric.normal_vector)
        return self._sem_tensor

    @SEMTensor.setter
    def SEMTensor(self, value):
        self._sem_tensor = value

    @property
    def AccelerationVector(self):
        if self._acceleration_vector is None:
            self._calculate_kinematic_quantities()
        return self._acceleration_vector

    @AccelerationVector.setter
    def AccelerationVector(self, value):
        self._acceleration_vector = value

    @property
    def ExpansionScalar(self):
        if self._expansion_scalar is None:
            self._calculate_kinematic_quantities()
        return self._expansion_scalar

    @ExpansionScalar.setter
    def ExpansionScalar(self, value):
        self._expansion_scalar = value

    @property
    def ShearTensor(self):
        if self._shear_tensor is None:
            self._calculate_kinematic_quantities()
        return self._shear_tensor

    @ShearTensor.setter
    def ShearTensor(self, value):
        self._shear_tensor = value

    @property
    def VorticityTensor(self):
        if self._vorticity_tensor is None:
            self._calculate_kinematic_quantities()
        return self._vorticity_tensor

    @VorticityTensor.setter
    def VorticityTensor(self, value):
        self._vorticity_tensor = value

    @property
    def VorticityVector(self):
        if self._vorticity_vector is None:
            self._calculate_kinematic_quantities()
        return self._vorticity_vector

    @VorticityVector.setter
    def VorticityVector(self, value):
        self._vorticity_vector = value

    @property
    def GravitoElectricTensor(self):
        if self._E_tensor is None:
            self._E_tensor = GravitoElectricTensor.from_weyl(self.WeylTensor, self.NormalVector)
        return self._E_tensor

    @GravitoElectricTensor.setter
    def GravitoElectricTensor(self, value):
        self._E_tensor = value

    @property
    def GravitoMagneticTensor(self):
        if self._H_tensor is None:
            self._H_tensor = GravitoMagneticTensor.from_weyl(self.WeylTensor, self.NormalVector)
        return self._H_tensor

    @GravitoMagneticTensor.setter
    def GravitoMagneticTensor(self, value):
        self._H_tensor = value


    def spatial_trace(self, S):
        """


        """
        h = self.ProjectorTensor.change_config("uu")
        r = tensor_product(h, S.change_config("ll"), 0, 0)
        r = tensorcontraction(r.arr, (0,1))
        return r


    def projection(self, T):
        """

        """
        h = self.ProjectorTensor
        for i in range(T.order):
            if T.config[i] == "l":
                h = h.change_config("lu")
            elif T.config[i] == "u":
                h = h.change_config("ll")
            T = tensor_product(h, T, 1, i)
        return OPTDecompositionTensor(T.arr, nvec=self.Metric.normal_vector, syms=T.syms, config=T.config, parent_metric=self.Metric)

    def projected_vector_dual(self, S):
        """

        """
        S = S.change_config("uu")
        S = S.antisymmetric_part()
        eps = self.ProjectedAlternatingTensor
        r = tensor_product(eps, S, 1, 0)
        r = tensorcontraction(r.arr, (1,2))/2
        return OPTDecompositionTensor(r, nvec=self.Metric.normal_vector, syms=S.syms, config="l", parent_metric=self.Metric)


    def pstf(self, S):
        """


        """
        h = self.ProjectorTensor
        S = S.change_config("ll")
        h = h.change_config("lu")
        l = tensor_product(h, S, 1, 0)
        l = tensor_product(h, l, 1, 1)
        l = l.symmetric_part().arr
        r = tensor_product(h.change_config("uu"), S, 0, 0)
        r = tensor_product(r, h.change_config("ll"))
        r = tensorcontraction(r.arr, (0,1))
        return OPTDecompositionTensor( l - r/3, nvec=self.Metric.normal_vector, syms=S.syms, config="ll", parent_metric=self.Metric)


    def vector_product(self, V, W):
        """

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
                r =  r.symmetric_part(0,1)
                r = r.arr
        elif V.order == 2:
            if W.order == 1:
                return self.vector_product(W, V)
            S = V.change_config("ul")
            Q = W.change_config("uu")
            r = tensor_product(eps, S, 1, 0)
            r = tensor_product(r, Q, 2, 0)
            r = tensorcontraction(r.arr, (1, 2))
        else:
            raise Exception("Higher orders not implemented")
        return OPTDecompositionTensor( r, nvec=V.normal_vector, syms=V.syms,
                                        config="ll" if V.order == 1 and W.order == 2 else "l", parent_metric=V.parent_metric)



    def projected_covariant_derivative(self, T):
        """


        """
        NablaT = self.covariant_derivative(T)
        u = self.Metric.normal_vector.change_config("u")
        h = self.ProjectorTensor
        Tdot = tensor_product(u, NablaT, 0, 0)
        DT = tensor_product(h.change_config("lu"), NablaT, 1, 0)
        for i in range(1, DT.order):
            if DT.config[i] == "l":
                DT = tensor_product(h.change_config("lu"), DT, 1, i)
            else:
                DT = tensor_product(h.change_config("ul"), DT, 1, i)
        return Tdot, DT


    def div(self, S):
        """

        """
        _, DS = self.projected_covariant_derivative(S)
        config = "u" + DS.config[1:]
        DS = DS.change_config(config)
        if S.order == 1:
            return tensorcontraction(DS.tensor(), (0,1))
        else:
            return OPTDecompositionTensor( tensorcontraction(DS.tensor(), (0, S.order)),
                                            nvec=S.normal_vector, syms=S.syms, config=config[1:-1], parent_metric=self.Metric)

    def curl(self, S):
        """

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
            return OPTDecompositionTensor( tensorcontraction( tensor_product( self.ProjectedAlternatingTensor, DS, 1, 0).arr, (1,2)),
                                            nvec=S.normal_vector, syms=S.syms, config="l", parent_metric=self.Metric)
        elif S.order == 2:
            return OPTDecompositionTensor( tensorcontraction( tensor_product( self.ProjectedAlternatingTensor, DS, 0, 0).arr, (0,3)),
                                            nvec=S.normal_vector, syms=S.syms, config="ll", parent_metric=self.Metric).symmetric_part()

    def temporal_rotation(self, S):
        """

        """
        u = self.Metric.normal_vector.change_config("u")
        eps_dot = self.ProjectedAlternatingTensorDot
        if S.order == 1:
            r = tensor_product(eps_dot, S.change_config("u"), 1, 0)
            r = tensor_product(u, r, 0, 1)
            return OPTDecompositionTensor(-r.arr, nvec=u, syms=S.syms, config="l", parent_metric=self.Metric)
        elif S.order == 2:
            r = tensor_product(eps_dot, S.change_config("lu"), 1, 1)
            r = tensor_product(u, r, 0, 0).symmetric_part()
            return OPTDecompositionTensor(-r.arr, nvec=u, syms=S.syms, config="ll", parent_metric=self.Metric)
        else:
            raise ValueError("Higher Orders not implemented")

