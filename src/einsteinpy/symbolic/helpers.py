import warnings

import numpy as np
import sympy
from sympy import ImmutableDenseNDimArray, derive_by_array


def raise_warning(WarningType, message):
    warnings.warn(message, WarningType)


def _flatten_list(seq):
    # flatten an arbitarily nested list
    if not seq:
        return []
    if not isinstance(seq[0], list):
        return [seq[0]] + _flatten_list(seq[1:])
    return _flatten_list(seq[0]) + _flatten_list(seq[1:])


def simplify_sympy_array(arr):
    """
    Function to simplify sympy expression or array.

    This function is explicitely defined as native ``simplify`` function within sympy
    stopped working with sympy version change.

    Parameters
    ----------
    arr : ~sympy.tensor.array.ndim_array.NDimArray or ~sympy.core.expr.Expr
        Any sympy array or expression.

    Returns
    -------
    sympy.tensor.array.ndim_array.NDimArray or ~sympy.core.expr.Expr
        Simplified sympy array or expression.

    """
    try:
        flattened_list = _flatten_list(arr.tolist())
        simplified_flattened_list = [sympy.simplify(e) for e in flattened_list]
        return sympy.Array(simplified_flattened_list, arr.shape)
    except AttributeError:
        return sympy.simplify(arr)
    except IndexError:
        return sympy.Array(sympy.simplify(sum(arr)))


def sympy_to_np_array(arr):
    """
    Function to convert sympy to numpy array

    Parameters
    ----------
    arr : ~sympy.tensor.array.ndim_array.NDimArray
        Sympy Array

    Returns
    -------
    ~numpy.ndarray
        Numpy Array

    """
    return np.array(arr.tolist()).reshape(arr.shape)


class TransformationMatrix(ImmutableDenseNDimArray):
    """
    Class for defining transformation matrix for basis change of vectors and tensors.

    """

    def __init__(
        self,
        iterable,
        old_coords,
        new_coords,
        shape=None,
        old2new=None,
        new2old=None,
        **kwargs
    ):
        """
        Constructor.

        Parameters
        ----------
        iterable : iterable-object
            2D list or array to pass a matrix.
        old_coords : list or tuple
            list of old coordinates.
            For example, ``[x, y]``.
        new_coords : list or tuple
            list of new coordinates.
            For example, ``[r, theta]``.
        shape : tuple, optional
            shape of the transformation matrix. Usually, not required.
            Defaults to ``None``.
        old2new : list or tuple, optional
            List of expressions for new coordinates in terms of old coordinates.
            For example, ``[x**2+y**2, atan2(y, x)]``.
        new2old : list or tuple, optional
            List of expressions for old coordinates in terms of new coordinates.
            For example, ``[r*cos(theta), r*sin(theta)]``.

        Raises
        ------
        ValueError
            Raised when tensor has a rank not equal to 2.
            This is because, a matrix is expected.

        """

        # __new__() is called after __init__() automatically
        self._inv = None
        self.old_coords, self.new_coords = old_coords, new_coords
        self.old2new, self.new2old = old2new, new2old

    def __new__(
        cls,
        iterable,
        old_coords,
        new_coords,
        shape=None,
        old2new=None,
        new2old=None,
        **kwargs
    ):
        obj = super(TransformationMatrix, cls).__new__(cls, iterable, shape, **kwargs)

        if not obj.rank() == 2:
            raise ValueError("Expected a tensor with rank 2")

        return obj

    @classmethod
    def from_new2old(cls, old_coords, new_coords, new2old, **kwargs):
        """
        Classmethod to obtain transformation matrix from old coordinates expressed
        as a function of new coordinates.

        Parameters
        ----------
        old_coords : list or tuple
            list of old coordinates.
            For example, ``[x, y]``.
        new_coords : list or tuple
            list of new coordinates.
            For example, ``[r, theta]``.
        new2old : list or tuple, optional
            List of expressions for old coordinates in terms of new coordinates.
            For example, ``[r*cos(theta), r*sin(theta)]``.

        """
        tmp_array = derive_by_array(new2old, new_coords)
        tmp_array = sympy_to_np_array(tmp_array)
        derivative_array = np.reciprocal(tmp_array)
        return cls(
            derivative_array,
            old_coords,
            new_coords,
            old2new=None,
            new2old=new2old,
            **kwargs
        )

    @classmethod
    def from_old2new(cls, old_coords, new_coords, old2new=None, new2old=None, **kwrags):
        raise NotImplementedError

    def inv(self):
        """
        Returns inverse of the transformation matrix

        Returns
        -------
        ~sympy.tensor.array.ndim_array.NDimArray
            Inverse of the matrix.

        """
        if self._inv is not None:
            return self._inv

        tmp_array = sympy_to_np_array(self)
        tmp_array = np.reciprocal(tmp_array)
        self._inv = ImmutableDenseNDimArray(tmp_array)
        return self._inv


def _change_name(curr_name: str, context: str) -> str:
    """
    Function to add descriptive tags to Tensor name, after the tensor is modified.
    Currently handles Lorentz Transformation and Config Changes.

    Parameters
    ----------
    curr_name : str
            Current name of the tensor
    context : str
            Context of name change - '__lt', for lorentz_transformation
                                   - '__' + newconfig, for config_change (cc)
    Returns
    -------
    str
        Altered name of the tensor, with appropriate tags

    """
    return curr_name + context if (curr_name is not None) else None
