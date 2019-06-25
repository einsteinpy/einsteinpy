import numpy as np
import sympy

from einsteinpy.symbolic import vaccum_metrics, MetricTensor


def schwarzschild_metric():
    symbolstr = "t r theta phi"
    syms = sympy.symbols(symbolstr)
    sch = vaccum_metrics.SchwarzschildMetric(syms)
    return sch


def kerr_metric():
    symbolstr = "t r theta phi"
    syms = sympy.symbols(symbolstr)
    sch = vaccum_metrics.KerrMetric(syms)
    return sch

