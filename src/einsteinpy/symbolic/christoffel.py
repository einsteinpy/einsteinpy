import numpy as np
import sympy

from einsteinpy.symbolic.helpers import _change_name
from einsteinpy.symbolic.tensor import BaseRelativityTensor, _change_config


class ChristoffelSymbols(BaseRelativityTensor):
    """
    Class for defining Christoffel symbols.

    .. note::
        Christoffel symbols are NOT tensors. They are "connection coefficients"
        that describe the Levi-Civita connection of a given metric tensor.
        Unlike true tensors, Christoffel symbols do not transform as tensors
        under coordinate transformations - they transform in a more complex way
        that involves derivatives of the coordinate transformation.

        While this class inherits from ``BaseRelativityTensor`` for code reuse
        purposes, Christoffel symbols are fundamentally different from tensors.
        The mathematical reason is that tensors transform homogeneously under
        coordinate changes, while Christoffel symbols have an additional
        inhomogeneous term involving derivatives of the transformation.

        For this reason, you may want to use the ``matrix()`` method instead
        of ``tensor()`` when accessing the underlying array, to avoid the
        implication that Christoffel symbols are tensors.
    """

    def __init__(
        self, arr, syms, config="ull", parent_metric=None, name="ChristoffelSymbols"
    ):
        """
        Constructor and Initializer

        Parameters
        ----------
        arr : ~sympy.tensor.array.dense_ndim_array.ImmutableDenseNDimArray or list
            Sympy Array or multi-dimensional list containing Sympy Expressions
        syms : tuple or list
            Tuple of crucial symbols denoting time-axis, 1st, 2nd, and 3rd axis (t,x1,x2,x3)
        config : str
            Configuration of contravariant and covariant indices in tensor. 'u' for upper and 'l' for lower indices. Defaults to 'ull'.
        parent_metric : ~einsteinpy.symbolic.metric.MetricTensor
            Metric Tensor from which Christoffel symbol is calculated. Defaults to None.
        name : str
            Name of the Christoffel Symbols Tensor. Defaults to "ChristoffelSymbols".

        Raises
        ------
        TypeError
            Raised when arr is not a list or sympy Array
        TypeError
            syms is not a list or tuple
        ValueError
            config has more or less than 3 indices

        """
        super(ChristoffelSymbols, self).__init__(
            arr=arr, syms=syms, config=config, parent_metric=parent_metric, name=name
        )
        self._order = 3
        if not len(self.config) == self._order:
            raise ValueError("config should be of length {}".format(self._order))

    @classmethod
    def from_metric(cls, metric):
        """
        Get Christoffel symbols calculated from a metric tensor

        Parameters
        ----------
        metric : ~einsteinpy.symbolic.metric.MetricTensor
            Space-time Metric from which Christoffel Symbols are to be calculated

        """
        dims = metric.dims
        tmplist = np.zeros((dims, dims, dims), dtype=int).tolist()
        mat, syms = metric.lower_config().tensor(), metric.symbols()
        matinv = sympy.Matrix(mat.tolist()).inv()
        for t in range(dims**3):
            # i,j,k each goes from 0 to (dims-1)
            # hack for codeclimate. Could be done with 3 nested for loops
            k = t % dims
            j = (int(t / dims)) % (dims)
            i = (int(t / (dims**2))) % (dims)
            if k <= j:
                tmpvar = 0
                for n in range(dims):
                    tmpvar += (matinv[i, n] / 2) * (
                        sympy.diff(mat[n, j], syms[k])
                        + sympy.diff(mat[n, k], syms[j])
                        - sympy.diff(mat[j, k], syms[n])
                    )
                tmplist[i][j][k] = tmplist[i][k][j] = tmpvar
        return cls(tmplist, syms, config="ull", parent_metric=metric)

    def change_config(self, newconfig="lll", metric=None):
        """
        Changes the index configuration(contravariant/covariant)

        Parameters
        ----------
        newconfig : str
            Specify the new configuration. Defaults to 'lll'
        metric : ~einsteinpy.symbolic.metric.MetricTensor or None
            Parent metric tensor for changing indices.
            Already assumes the value of the metric tensor from which it was initialized if passed with None.
            Compulsory if not initialized with 'from_metric'. Defaults to None.

        Returns
        -------
        ~einsteinpy.symbolic.christoffel.ChristoffelSymbols
            New tensor with new configuration. Defaults to 'lll'

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
        new_obj = ChristoffelSymbols(
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
            ~einsteinpy.symbolic.christoffel.ChristoffelSymbols
                lorentz transformed tensor

        """
        t = super(ChristoffelSymbols, self).lorentz_transform(transformation_matrix)
        return ChristoffelSymbols(
            t.tensor(),
            syms=self.syms,
            config=self._config,
            parent_metric=None,
            name=_change_name(self.name, context="__lt"),
        )

    def matrix(self):
        """
        Returns the Christoffel symbols as a sympy Array.

        This method is provided as an alternative to ``tensor()`` to emphasize
        that Christoffel symbols are not tensors but rather connection coefficients.

        Returns
        -------
        ~sympy.tensor.array.dense_ndim_array.ImmutableDenseNDimArray
            Sympy Array object representing the Christoffel symbols

        """
        return self.arr
