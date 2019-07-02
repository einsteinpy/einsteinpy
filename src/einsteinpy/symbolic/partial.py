from sympy import Basic, Expr, Array, Mul, S, Symbol, diff
from sympy.tensor.tensor import Tensor as SympyTensor
from sympy.core.decorators import call_highest_priority
from .tensor import Tensor


class DiffOperator(Expr):
    _op_priority = 11
    is_commutative = False

    def __new__(cls, *args, left=None):
        if left == S.Zero:
            return S.Zero
        obj = Basic.__new__(cls, *args)
        obj.left = left if left is not None else S.One
        return obj

    def __repr__(self):
        syms = [s.name for s in self.args]
        name = "\u2202(%s)" % ",".join(syms)
        return str(Mul(self.left, Symbol(name, commutative=False)))

    def __str__(self):
        return self.__repr__()

    @call_highest_priority("__rmul__")
    def __mul__(self, other):
        if isinstance(other, DiffOperator):
            args = self.args + other.args
            return DiffOperator(*args)
        return Mul(self.left, diff(other, *self.args)).doit()

    @call_highest_priority("__mul__")
    def __rmul__(self, other):
        return DiffOperator(*self.args, left=other)


class PartialDerivative(Tensor):
    is_TensorDerivative = True

    def __new__(cls, metric, **kwargs):
        basis = list(map(DiffOperator, metric.coords))
        return super().__new__(cls, "\u2202", basis, metric, covar=(-1,))

    def __repr__(self):
        return self._print()

    def __str__(self):
        return self.__repr__()
