import numpy as np
from collections import namedtuple
from sympy import Array, Symbol, eye, ones, diag, symbols, simplify
from sympy.tensor.tensor import TensorIndexType, TensorIndex, TensorHead, \
    TensorType, TensExpr, tensorsymmetry
from sympy.tensor.tensor import Tensor as SympyTensor
from sympy.tensor.array import tensorcontraction, tensorproduct, permutedims
from sympy.core.compatibility import string_types

class AbstractTensor(object):
    """
    Wrapper class for sympy Array with attributes used for identification.
    """

    is_Tensor = True
    is_Metric = False
    is_Spacetime = False
    _array = None
    _inverse = None

    def __new__(cls, obj, matrix):
        obj._array = matrix
        return obj

    def as_matrix(self):
        """
        Return the data stored in the tensor as an instance of ~sympy.Matrix.

        Notes
        -----
        The ``tomatrix`` attribute required by this method will throw an error
        if the tensor is not of rank 2.

        """
        return self._array.tomatrix()

    def as_array(self):
        """
        Return the data stored in the tensor as an instance of ~sympy.Array.
        """
        return self._array.copy()

    def as_inverse(self):
        if self._inverse == None:
            self._inverse = Array(self.as_matrix().inv())
        return self._inverse

class IndexedTensor(AbstractTensor, SympyTensor):
    """
    Class representing a Tensor that has been evaluated with indices.

    Used as the base object for creating algebraic expressions with tensors.
    Inherited are methods such as __mul__ to allow sympy to manage dummy
    indices when multiplied with other tensors.

    Generated when a Tensor is called as a function with indices as arguments.

    """

    def __new__(cls, tensor, indices, **kwargs):
        obj = SympyTensor.__new__(cls, tensor, indices, **kwargs)
        metrics = map(getattr, indices, len(indices)*['tensor_index_type'])
        repl = {}
        repl[obj] = tensor.covariance_transform(*indices)
        for metric in metrics:
            repl[metric] = metric.as_array()
        obj.replacement_dict = repl
        return AbstractTensor.__new__(cls, obj, tensor._array)

class Tensor(AbstractTensor, TensorHead):
    """
    Class for representing a Sympy TensorHead object that have an associated
    array of data elements/expressions to be substituted when requested.

    """

    def __new__(cls, symbol, matrix, metric, **kwargs):
        """
        Create a new Tensor object.

        Parameters
        ----------
        symbol : str
            Name of the tensor and the symbol to denote it by when printed.
        matrix : (list, tuple, ~sympy.Matrix, ~sympy.Array)
            Matrix representation of the tensor to be used in substitution.
            Can be of any type that is acceptable by ~sympy.Array.
        metric : Metric
            Classify the tensor as being defined in terms of a metric.

        Notes
        -----
        If the parameter ``symmetry`` is passed, the tensor object will defined
        using a specific symmetry. Example values are (see sympy documentation
        for the function ``tensorsymmetry``):
        ``[[1]]``         vector
        ``[[1]*n]``       symmetric tensor of rank ``n``
        ``[[n]]``         antisymmetric tensor of rank ``n``
        ``[[2, 2]]``      monoterm slot symmetry of the Riemann tensor
        ``[[1],[1]]``     vector*vector
        ``[[2],[1],[1]]`` (antisymmetric tensor)*vector*vector

        Examples
        --------
        >>> from sympy import diag, symbols
        >>> from einsteinpy.symbolic.tensor import *
        >>> E1, E2, E3, B1, B2, B3 = symbols('E1:4 B1:4')
        >>> em = [[0, -E1, -E2, -E3],
                  [E1, 0, -B3, B2],
                  [E2, B3, 0, -B1],
                  [E3, -B2, B1, 0]]
        >>> eta = SpacetimeMetric('eta', diag(1, -1, -1, -1))
        >>> F = Tensor('F', em, eta, symmetry=[[2]])
        >>> mu, nu = indices('mu nu', eta)
        >>> expr = F(mu, nu) + F(nu, mu)
        >>> expand_tensor(expr)
        0
        >>> expr = F(mu, nu) * F(-mu, -nu)
        >>> expand_tensor(expr)
        2*B_1**2 + 2*B_2**2 + 2*B_3**2 - 2*E_1**2 - 2*E_2**2 - 2*E_3**2

        """
        array = Array(matrix)
        sym = kwargs.pop('symmetry', [[1]*array.rank()])
        sym = tensorsymmetry(*sym)
        symtype = TensorType(array.rank()*[metric], sym)
        obj = TensorHead.__new__(cls, symbol, symtype, **kwargs)
        # resolves a bug with pretty printing.
        # TODO: Consider renaming this class to avoid conflicts.
        obj.__class__.__name__ = 'TensorHead'
        return AbstractTensor.__new__(cls, obj, array)

    def __repr__(self):
        return self._print()

    def __str__(self):
        return str(self.__repr__())

    def __call__(self, *args, **kwargs):
        new = IndexedTensor(self, args, **kwargs)
        return new.doit()

    def covariance_transform(self, *indices):
        """
        Return the array associated with this tensor with indices set according
        to arguments.

        Parameters
        ----------
        indices : TensorIndex
            Defines the covariance and contravariance of the returned array.

        Examples
        --------
        >>> from sympy import diag, symbols, sin
        >>> from einsteinpy.symbolic.tensor import *
        >>> r, th = symbols('r theta', real=True)
        >>> schwarzschild = diag(1-1/r, -1/(1-1/r), -r**2, -r**2*sin(th)**2)
        >>> g = SpacetimeMetric('g', schwarzschild)
        >>> mu, nu = indices('mu nu', g)
        >>> g.covariance_transform(-mu, -nu)
        [[1/(1 - 1/r), 0, 0, 0], [0, 1 - 1/r, 0, 0], [0, 0, -1/r**2, 0], [0, 0, 0, -1/(r**2*sin(theta)**2)]]

        """
        array = self.as_array()
        for pos, idx in enumerate(indices):
            if not idx.is_up:
                metric = idx.tensor_index_type.metric.as_inverse()
                new = tensorcontraction(tensorproduct(metric, array), (1, 2+pos))
                permu = list(range(len(indices)))
                permu[0], permu[pos] = permu[pos], permu[0]
                array = permutedims(new, permu)
        return array

class Metric(AbstractTensor, TensorIndexType):
    """
    Class representing a tensor that raises and lowers indices.
    """

    _MetricId = namedtuple('MetricId', ['name', 'antisym'])
    is_Metric = True

    def __new__(cls, symbol, matrix, **kwargs):
        """
        Create a new Metric object.

        Parameters
        ----------
        symbol : str
            Name of the tensor and the symbol to denote it by when printed.
        matrix : (list, tuple, ~sympy.Matrix, ~sympy.Array)
            Matrix representation of the tensor to be used in substitution.
            Can be of any type that is acceptable by ~sympy.Array.

        Examples
        --------
        >>> from sympy import diag
        >>> from einsteinpy.symbolic.tensor import *
        >>> eta = SpacetimeMetric('eta', diag(1, -1, -1, -1))
        >>> mu, nu = indices('mu nu', eta)
        >>> expr = eta(mu, nu) * eta(-mu, -nu)
        >>> expand_tensor(expr)
        4

        """

        array = Array(matrix)
        if array.rank() != 2 or array.shape[0] != array.shape[1]:
            raise ValueError(f'matrix must be square, received matrix of shape {array.shape}')
        obj = TensorIndexType.__new__(cls, symbol,
                                      metric=cls._MetricId(symbol, False),
                                      dim=array.shape[0], dummy_fmt=symbol,
                                      **kwargs)
        obj.metric = Tensor(obj.name, array, obj)
        return AbstractTensor.__new__(cls, obj, array)

    def __getattr__(self, attr):
        if hasattr(self.metric, attr):
            return getattr(self.metric, attr)
        return TensorIndexType.__getattribute__(self, attr)

    def __call__(self, *args):
        return self.metric(*args)

class SpacetimeMetric(Metric):
    """
    Class representing psuedo-Riemannian metrics.
    """

    is_Spacetime = True

    def __new__(cls, symbol, matrix, timelike=True, **kwargs):
        obj = super().__new__(cls, symbol, matrix, **kwargs)
        if obj.dim > 4:
            raise ValueError('metrics on spacetime must be at most 4-dimensional')
        obj.is_timelike = timelike
        obj.is_spacelike = not timelike
        return obj

    def reverse_signature(self):
        self._array *= -1
        self._replacement_dict = {self : self._array}
        self.is_timelike = not self.is_timelike
        self.is_spacelike = not self.is_spacelike
        return self.signature

    @property
    def signature(self):
        sign = -1 if self.is_timelike else 1
        sig = sign * ones(1, self.dim)
        sig[0] *= -1
        return tuple(sig)

class Index(TensorIndex):
    """
    Class for a symbolic representation of a tensor index with respect to a metric.
    """

    def __new__(cls, symbol, metric, is_up=True, **kwargs):
        return super().__new__(cls, symbol, metric, is_up=is_up, **kwargs)

    def __neg__(self):
        return Index(self.name, self.tensor_index_type, (not self.is_up))

def indices(s, metric=Euclidean(), is_up=True):
    if isinstance(s, string_types):
        a = [x.name for x in symbols(s, seq=True)]
    else:
        raise ValueError(f'expected a string, received object of type {type(s)}')
    idxs = [Index(idx, metric, is_up) for idx in a]
    if len(idxs) == 1:
        return idxs[0]
    return idxs

def get_replacement_dict(expr):
    repl = {}
    terms = expr.args if not isinstance(expr, IndexedTensor) else [expr]
    for term in terms:
        repl.update(term.replacement_dict)
    return repl

def expand_tensor(expr, idxs=None):
    """
    Evaluate a tensor expression and return the resultant array.
    """

    repl = get_replacement_dict(expr)
    for m in filter(lambda m: m.is_Metric, repl.keys()):
        expr = expr.contract_metric(m.metric)
        if not isinstance(expr, TensExpr):
            return expr
    if idxs is None:
        idxs = expr.get_free_indices()
    return expr.replace_with_arrays(repl, idxs)
