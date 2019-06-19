import numpy as np
import sympy


class Tensor:
    """
    Base Class
    """

    def __init__(self, syms):
        """
        Constructor and Initializer
        :param syms:
        """
        self.syms = syms
        self.dims = len(self.syms)
        self.generic_list = (np.zeros(shape=(self.dims, self.dims, self.dims, self.dims), dtype=int)).tolist()

    def _simplify_tensor_helper(self, v):
        """

        :param v:
        :return:
        """
        self.returnval = None
        if isinstance(v, list):
            newlist = list()
            for t in v:
                newlist.append(self._simplify_tensor_helper(t))
            self.returnval = newlist
        else:
            self.returnval = sympy.simplify(v)
        return self.returnval

    def simplify_tensor(self, ndlist):
        """
        Returns a N-Dimensional list of simplified symbolic tensor

        Parameters
        ----------
        ndlist : list
            N-Dimensional list of sympy expressions representing some tensor

        Returns
        -------
        list
            N-Dimensional list

        """
        return self._simplify_tensor_helper(ndlist)
