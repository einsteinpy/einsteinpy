import itertools as it
from collections import defaultdict
from sympy import S, Array, symbols, diff, simplify
from sympy.tensor.tensor import TensorIndex, TensorHead, TensorType, \
    TensExpr, TensMul, TensAdd, TensorManager, tensorsymmetry
from sympy.tensor.tensor import Tensor as SympyTensor
from sympy.tensor.array import tensorcontraction, tensorproduct, permutedims
from sympy.core.compatibility import string_types

class _ReplacementManager(dict):
    def has(self, tensor):
        for key in self.keys():
            if key.args[0] == tensor.args[0]:
                return True
        return False

    def remove(self, tensor):
        for key in self.keys():
            if key.args[0] == tensor.args[0]:
                self.pop(key)
                break

    def __setitem__(self, tensor, array):
        if not self.has(tensor):
            self.update({tensor: array})

ReplacementManager = _ReplacementManager()

class AbstractTensor(object):
    """
    Wrapper class for sympy Array with attributes used for identification.
    """
    is_Tensor = True
    is_Metric = False
    is_Spacetime = False
    is_TensorDerivative = False
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
        array = tensor.covariance_transform(*indices)
        ReplacementManager[obj] = array
        metrics = [idx.tensor_index_type for idx in indices]
        for metric in metrics:
            ReplacementManager[metric] = metric.as_array()
        return AbstractTensor.__new__(cls, obj, array)

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
        >>> t, x, y, z = symbols('t x y z')
        >>> eta = Metric('eta', [t, x, y, z], diag(1, -1, -1, -1))
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
        comm = kwargs.pop('comm', 'general')
        covar = tuple(kwargs.pop('covar', array.rank()*[1]))
        if len(covar) != array.rank():
            raise ValueError(f'covariance signature {covar} does not match tensor rank {array.rank()}')
        obj = TensorHead.__new__(cls, symbol, symtype, comm=comm, **kwargs)
        # resolves a bug with pretty printing.
        # TODO: Consider renaming this class to avoid conflicts.
        obj.__class__.__name__ = 'TensorHead'
        obj.covar = covar
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
        >>> t, r, th, ph = symbols('t r theta phi')
        >>> schwarzschild = diag(1-1/r, -1/(1-1/r), -r**2, -r**2*sin(th)**2)
        >>> g = Metric('g', [t, r, th, ph], schwarzschild)
        >>> mu, nu = indices('mu nu', g)
        >>> g.covariance_transform(mu, nu)
        [[1/(1 - 1/r), 0, 0, 0], [0, 1 - 1/r, 0, 0], [0, 0, -1/r**2, 0], [0, 0, 0, -1/(r**2*sin(theta)**2)]]

        """
        array = self.as_array()
        for pos, idx in enumerate(indices):
            if idx.is_up ^ (self.covar[pos] > 0):
                if idx.is_up:
                    metric = idx.tensor_index_type.metric.as_inverse()
                else:
                    metric = idx.tensor_index_type.metric.as_array()
                new = tensorcontraction(tensorproduct(metric, array), (1, 2 + pos))
                permu = list(range(len(indices)))
                permu[0], permu[pos] = permu[pos], permu[0]
                array = permutedims(new, permu)
        return array

    def simplify(self):
        array = simplify(self.as_array())
        self._array = array
        ReplacementManager.remove(self)
        return array

class Index(TensorIndex):
    """
    Class for a symbolic representation of a tensor index with respect to a metric.
    """
    def __new__(cls, symbol, metric, is_up=True, **kwargs):
        return super().__new__(cls, symbol, metric, is_up=is_up, **kwargs)

    def __neg__(self):
        return Index(self.name, self.tensor_index_type, (not self.is_up))

def expand_tensor(expr, idxs=None):
    """
    Evaluate a tensor expression and return the resulting array.
    """
    if idxs is None:
        idxs = TensMul(expr).get_free_indices()
    return expr.replace_with_arrays(ReplacementManager, idxs)

def indices(s, metric, is_up=True):
    if isinstance(s, string_types):
        a = [x.name for x in symbols(s, seq=True)]
    else:
        raise ValueError(f'expected a string, received object of type {type(s)}')
    idxs = [Index(idx, metric, is_up) for idx in a]
    if len(idxs) == 1:
        return idxs[0]
    return idxs

TensorManager.set_comm('general', 'metric', 0)
TensorManager.set_comm('metric', 'metric', 0)
TensorManager.set_comm('general', 'general', 0)
TensorManager.set_comm('partial', 'partial', 0)
