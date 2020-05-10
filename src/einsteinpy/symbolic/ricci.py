import sympy
from sympy import tensorcontraction, tensorproduct

from einsteinpy.symbolic.christoffel import ChristoffelSymbols
from einsteinpy.symbolic.helpers import _change_name, simplify_sympy_array
from einsteinpy.symbolic.metric import MetricTensor
from einsteinpy.symbolic.riemann import RiemannCurvatureTensor
from einsteinpy.symbolic.tensor import BaseRelativityTensor, _change_config


class RicciTensor(BaseRelativityTensor):
    """
    Class for defining Ricci Tensor
    """

    def __init__(self, arr, syms, config="ll", parent_metric=None, name="RicciTensor"):
        """
        Constructor and Initializer

        Parameters
        ----------
        arr : ~sympy.tensor.array.dense_ndim_array.ImmutableDenseNDimArray or list
            Sympy Array or multi-dimensional list containing Sympy Expressions
        syms : tuple or list
            Tuple of crucial symbols denoting time-axis, 1st, 2nd, and 3rd axis (t,x1,x2,x3)
        config : str
            Configuration of contravariant and covariant indices in tensor. 'u' for upper and 'l' for lower indices. Defaults to 'll'.
        parent_metric : ~einsteinpy.symbolic.metric.MetricTensor or None
            Corresponding Metric for the Ricci Tensor.
            Defaults to None.
        name : str
            Name of the Tensor. Defaults to "RicciTensor".

        Raises
        ------
        TypeError
            Raised when arr is not a list or sympy Array
        TypeError
            syms is not a list or tuple
        ValueError
            config has more or less than 2 indices

        """
        super(RicciTensor, self).__init__(
            arr=arr, syms=syms, config=config, parent_metric=parent_metric, name=name
        )
        self._order = 2
        if not len(config) == self._order:
            raise ValueError("config should be of length {}".format(self._order))

    @classmethod
    def from_riemann(cls, riemann, parent_metric=None):
        """
        Get Ricci Tensor calculated from Riemann Tensor

        Parameters
        ----------
        riemann : ~einsteinpy.symbolic.riemann.RiemannCurvatureTensor
           Riemann Tensor
        parent_metric : ~einsteinpy.symbolic.metric.MetricTensor or None
            Corresponding Metric for the Ricci Tensor.
            None if it should inherit the Parent Metric of Riemann Tensor.
            Defaults to None.

        """
        if not riemann.config == "ulll":
            riemann = riemann.change_config(newconfig="ulll", metric=parent_metric)
        if parent_metric is None:
            parent_metric = riemann.parent_metric
        return cls(
            simplify_sympy_array(sympy.tensorcontraction(riemann.tensor(), (0, 2))),
            riemann.syms,
            config="ll",
            parent_metric=parent_metric,
        )

    @classmethod
    def from_christoffels(cls, chris, parent_metric=None):
        """
        Get Ricci Tensor calculated from Christoffel Tensor

        Parameters
        ----------
        chris : ~einsteinpy.symbolic.christoffel.ChristoffelSymbols
            Christoffel Tensor
        parent_metric : ~einsteinpy.symbolic.metric.MetricTensor or None
            Corresponding Metric for the Ricci Tensor.
            None if it should inherit the Parent Metric of Christoffel Symbols.
            Defaults to None.

        """
        rt = RiemannCurvatureTensor.from_christoffels(
            chris, parent_metric=parent_metric
        )
        return cls.from_riemann(rt)

    @classmethod
    def from_metric(cls, metric):
        """
        Get Ricci Tensor calculated from Metric Tensor

        Parameters
        ----------
        metric : ~einsteinpy.symbolic.metric.MetricTensor
            Metric Tensor

        """
        ch = ChristoffelSymbols.from_metric(metric)
        return cls.from_christoffels(ch, parent_metric=None)

    def change_config(self, newconfig="ul", metric=None):
        """
        Changes the index configuration(contravariant/covariant)

        Parameters
        ----------
        newconfig : str
            Specify the new configuration. Defaults to 'ul'
        metric : ~einsteinpy.symbolic.metric.MetricTensor or None
            Parent metric tensor for changing indices.
            Already assumes the value of the metric tensor from which it was initialized if passed with None.
            Compulsory if somehow does not have a parent metric. Defaults to None.

        Returns
        -------
        ~einsteinpy.symbolic.ricci.RicciTensor
            New tensor with new configuration. Defaults to 'ul'

        Raises
        ------
        Exception
            Raised when a parent metric could not be found.

        """
        if metric is None:
            metric = self._parent_metric
        if metric is None:
            raise Exception("Parent Metric not found, can't do configuration change")
        new_tensor = _change_config(self, metric, newconfig)
        new_obj = RicciTensor(
            new_tensor,
            self.syms,
            config=newconfig,
            parent_metric=metric,
            name=_change_name(self.name, context="__" + newconfig),
        )
        return new_obj

    def lorentz_transform(self, transformation_matrix):
        """
        Performs a Lorentz transform on the tensor.

        Parameters
        ----------
            transformation_matrix : ~sympy.tensor.array.dense_ndim_array.ImmutableDenseNDimArray or list
                Sympy Array or multi-dimensional list containing Sympy Expressions

        Returns
        -------
            ~einsteinpy.symbolic.ricci.RicciTensor
                lorentz transformed tensor

        """
        t = super(RicciTensor, self).lorentz_transform(transformation_matrix)
        return RicciTensor(
            t.tensor(),
            syms=self.syms,
            config=self._config,
            parent_metric=None,
            name=_change_name(self.name, context="__lt"),
        )


class RicciScalar(BaseRelativityTensor):
    """
    Class for defining Ricci Scalar
    """

    def __init__(self, arr, syms, parent_metric=None):
        """
        Constructor and Initializer

        Parameters
        ----------
        arr : ~sympy.tensor.array.dense_ndim_array.ImmutableDenseNDimArray or list
            Sympy Array, multi-dimensional list containing Sympy Expressions, or Sympy Expressions or int or float scalar
        syms : tuple or list
            Tuple of crucial symbols denoting time-axis, 1st, 2nd, and 3rd axis (t,x1,x2,x3)
        parent_metric : ~einsteinpy.symbolic.metric.MetricTensor or None
            Corresponding Metric for the Ricci Scalar.
            Defaults to None.

        Raises
        ------
        TypeError
            Raised when syms is not a list or tuple

        """
        super(RicciScalar, self).__init__(
            arr=arr,
            syms=syms,
            config="",
            parent_metric=parent_metric,
            name="RicciScalar",
        )
        self._order = 0

    @property
    def expr(self):
        """
        Retuns the symbolic expression of the Ricci Scalar
        """
        val = sum(self.arr)  # sympy not allowing indexing, temporary fix!
        return val

    @classmethod
    def from_riccitensor(cls, riccitensor, parent_metric=None):
        """
        Get Ricci Scalar calculated from Ricci Tensor

        Parameters
        ----------
        riccitensor: ~einsteinpy.symbolic.metric.RicciTensor
            Ricci Tensor
        parent_metric : ~einsteinpy.symbolic.metric.MetricTensor or None
            Corresponding Metric for the Ricci Scalar.
            Defaults to None.

        """

        if not riccitensor.config == "ul":
            riccitensor = riccitensor.change_config(
                newconfig="ul", metric=parent_metric
            )
        if parent_metric is None:
            parent_metric = riccitensor.parent_metric
        ricci_scalar = tensorcontraction(riccitensor.tensor(), (0, 1))
        return cls(
            simplify_sympy_array(ricci_scalar),
            riccitensor.syms,
            parent_metric=parent_metric,
        )

    @classmethod
    def from_riemann(cls, riemann, parent_metric=None):
        """
        Get Ricci Scalar calculated from Riemann Tensor

        Parameters
        ----------
        riemann : ~einsteinpy.symbolic.riemann.RiemannCurvatureTensor
           Riemann Tensor
        parent_metric : ~einsteinpy.symbolic.metric.MetricTensor or None
            Corresponding Metric for the Ricci Scalar.
            Defaults to None.

        """

        cg = RicciTensor.from_riemann(riemann, parent_metric=parent_metric)
        return cls.from_riccitensor(cg)

    @classmethod
    def from_christoffels(cls, chris, parent_metric=None):
        """
        Get Ricci Scalar calculated from Christoffel Tensor

        Parameters
        ----------
        chris : ~einsteinpy.symbolic.christoffel.ChristoffelSymbols
            Christoffel Tensor
        parent_metric : ~einsteinpy.symbolic.metric.MetricTensor or None
            Corresponding Metric for the Ricci Scalar.
            Defaults to None.

        """
        rt = RiemannCurvatureTensor.from_christoffels(
            chris, parent_metric=parent_metric
        )
        return cls.from_riemann(rt)

    @classmethod
    def from_metric(cls, metric):
        """
        Get Ricci Scalar calculated from Metric Tensor

        Parameters
        ----------
        metric : ~einsteinpy.symbolic.metric.MetricTensor
            Metric Tensor

        """
        ch = ChristoffelSymbols.from_metric(metric)
        return cls.from_christoffels(ch, parent_metric=None)
