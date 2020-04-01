import numpy as np
import sympy

from einsteinpy.symbolic import ChristoffelSymbols, MetricTensor, vacuum_metrics


def test_SchwarzschildMetric():
    symbolstr = "t r theta phi"
    sch = vacuum_metrics.SchwarzschildMetric(symbolstr)
    chl = ChristoffelSymbols.from_metric(sch)
    mat = chl.tensor()
    syms = sympy.symbols(symbolstr)
    G, M, c, a = sympy.symbols("G M c a")
    assert (
        sympy.simplify(mat[2, 3, 3] - (-1 * sympy.cos(syms[2]) * sympy.sin(syms[2])))
        == 0
    )
    assert sympy.simplify(mat[3, 3, 1] - syms[1] / (syms[1] ** 2)) == 0
    assert (
        sympy.simplify(
            (mat[1, 1, 1].subs({a: (2 * G * M / (c ** 2))}))
            - (G * M / (2 * G * M * syms[1] - c ** 2 * syms[1] ** 2))
        )
        == 0
    )
    assert chl.symbols() == syms
