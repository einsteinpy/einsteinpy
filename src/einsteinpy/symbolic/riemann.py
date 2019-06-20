import sympy

from .christoffel import ChristoffelSymbols
from .tensor import Tensor


class RiemannCurvatureTensor(Tensor):
    def __init__(self, list2d, syms):
        """
        Constructor and Initializer
        :param list2d:
        :param syms:
        """
        super.__init__(syms)
        self.list2d = list2d
        dims = self.dims
        christs = ChristoffelSymbols(self.list2d, self.syms)
        self.riemann_list = self.create_riemannlist(
            dims, self.generic_list, christs.christlist
        )

    def create_riemannlist(self, dims, generic_list, christslist):
        _counterlist = [i for i in range(dims ** 4)]
        for i in _counterlist:
            # t,s,r,n each goes from 0 to (dims-1)
            # hack for codeclimate. Could be done with 4 nested for loops
            n = i % dims
            r = (int(i / dims)) % (dims)
            s = (int(i / (dims ** 2))) % (dims)
            t = (int(i / (dims ** 3))) % (dims)
            temp = sympy.diff(christslist[t][s][n], syms[r]) - sympy.diff(
                christslist[t][r][n], syms[s]
            )
            for p in range(dims):
                temp += (christslist[p][s][n] * christslist[t][p][r]) - christslist[p][
                    r
                ][n] * christslist[t][p][s]
            generic_list[t][s][r][n] = sympy.simplify(temp)
        return generic_list
