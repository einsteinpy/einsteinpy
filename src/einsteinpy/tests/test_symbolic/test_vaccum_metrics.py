import numpy as np
import sympy

from einsteinpy.symbolic import MetricTensor, vaccum_metrics


def schwarzschild_metric():
    symbolstr = "t r theta phi"
    syms = sympy.symbols(symbolstr)
    sch = vaccum_metrics.SchwarzschildMetric(syms)
    return sch
