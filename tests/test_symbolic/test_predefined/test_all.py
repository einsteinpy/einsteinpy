import numpy as np
import pytest
from sympy import Array, diff, log, simplify, sqrt, symbols, tanh

from einsteinpy.symbolic import (
    MetricTensor,
    RicciScalar,
    RicciTensor,
    constants,
    simplify_sympy_array,
)
from einsteinpy.symbolic.predefined import (
    AlcubierreWarp,
    AntiDeSitter,
    AntiDeSitterStatic,
    BarriolaVilekin,
    BertottiKasner,
    BesselGravitationalWave,
    CMetric,
    Davidson,
    DeSitter,
    Ernst,
    FriedmannLemaitreRobertsonWalker,
    Godel,
    JanisNewmanWinicour,
    Kerr,
    KerrNewman,
    KerrNewmanDeSitter,
    Minkowski,
    MinkowskiCartesian,
    MinkowskiPolar,
    ReissnerNordstorm,
    Schwarzschild,
)


@pytest.mark.parametrize(
    "metric_instance",
    [
        AlcubierreWarp(),
        AntiDeSitter(),
        AntiDeSitterStatic(),
        DeSitter(),
        Schwarzschild(),
        Schwarzschild(c=1, sch=2),
        Kerr(),
        KerrNewman(),
        Minkowski(),
        MinkowskiCartesian(),
        MinkowskiPolar(),
        CMetric(),
        Davidson(),
        Godel(),
        BesselGravitationalWave(),
        BarriolaVilekin(),
        BertottiKasner(),
        Ernst(),
        JanisNewmanWinicour(),
        ReissnerNordstorm(),
        FriedmannLemaitreRobertsonWalker(),
        KerrNewmanDeSitter(),
    ],
)
def test_all_predefined_metrics(metric_instance):
    assert isinstance(metric_instance, MetricTensor)


@pytest.mark.parametrize(
    "m1, m2",
    [
        (Schwarzschild(), Kerr(a=0)),  # Schwarzschild is a special case of Kerr
        (Kerr(), KerrNewman(Q=0)),  # Kerr is a special case of Kerr-Newman
        (
            ReissnerNordstorm(),
            KerrNewman(a=0),
        ),  # Reissner-Nordstorm is a special case of Kerr-Newman
        (
            Schwarzschild(),
            ReissnerNordstorm(Q=0),
        ),  # Schwarzschild is a special case of Reissner-Nordstorm
    ],
)
def test_check_two_metrics_are_equal(m1, m2):
    zero_arr = Array(np.zeros(shape=m1.tensor().shape, dtype=int))
    assert simplify_sympy_array(m1.tensor() - m2.tensor()) == zero_arr


def test_AlcubierreWarp_metric():
    # Regression test to ensure metric is same as in paper (Alcubierre 1994).
    x_s, sigma, R, v = symbols("x_s sigma R v")
    metric = AlcubierreWarp(x_s=x_s, sigma=sigma, R=R, v=v)
    _, x, y, z = metric.symbols()
    r_s = sqrt((x - x_s) ** 2 + y**2 + z**2)
    f = (tanh(sigma * (r_s + R)) - tanh(sigma * (r_s - R))) / (2 * tanh(sigma * R))
    expected = Array(
        [
            [v**2 * f**2 - 1, -v * f, 0, 0],
            [-v * f, 1, 0, 0],
            [0, 0, 1, 0],
            [0, 0, 0, 1],
        ]
    )
    assert metric.tensor() == expected


def test_DeSitter_has_constant_ricci_scalar():
    metric = DeSitter()
    al = symbols("alpha")
    R = RicciScalar.from_metric(metric).expr
    assert simplify(R - 12 / al**2) == 0


def test_Minkowski_equality():
    # Minkowski and MinkowskiCartesian are same
    assert simplify_sympy_array(MinkowskiCartesian().tensor()) == simplify_sympy_array(
        Minkowski().tensor()
    )


def test_JanisNewmanWinicour_scalar_field_and_Schwarzschild_limit():
    c, G = constants.c, constants.G
    gam, M, r = symbols("gam M r")
    # gam = 1 is Schwarzschild, which uses the opposite signature
    jnw = JanisNewmanWinicour(gam=1).tensor()
    sch = Schwarzschild(sch=2 * G * M / c**2).tensor()
    assert simplify_sympy_array(jnw + sch) == Array(np.zeros((4, 4), dtype=int))
    # A massless scalar field phi(r) only sources R_rr, which is
    # (1 - gam^2) / 2 * (d ln(alpha) / dr)^2
    ric = RicciTensor.from_metric(JanisNewmanWinicour()).tensor()
    alpha = 1 - 2 * G * M / (c**2 * gam * r)
    expected = (1 - gam**2) / 2 * diff(log(alpha), r) ** 2
    for i in range(4):
        for j in range(4):
            if (i, j) != (1, 1):
                assert simplify(ric[i, j]) == 0
    assert simplify(ric[1, 1] - expected) == 0


def test_BesselGravitationalWave_is_vacuum():
    # Kramer (1999) is a vacuum solution, so the Ricci tensor must vanish
    metric = BesselGravitationalWave()
    zero_arr = Array(np.zeros(shape=metric.tensor().shape, dtype=int))
    assert RicciTensor.from_metric(metric).tensor() == zero_arr
    # The wave depends on time through U = C J0(rho) cos(t)
    assert symbols("t") in metric.tensor()[3, 3].free_symbols
