"""
Physics checks for ``lorentz_transform`` across the symbolic module.

The main check builds the boosted metric by hand, recomputes each curvature
quantity from it in the new coordinates, and compares that with
``lorentz_transform`` of the quantity computed in the old frame.
"""

from functools import lru_cache

import pytest
import sympy
from sympy import ImmutableMatrix, Rational, cosh, simplify, sinh, symbols

from einsteinpy.symbolic import (
    ChristoffelSymbols,
    EinsteinTensor,
    GenericVector,
    MetricTensor,
    RicciScalar,
    RicciTensor,
    RiemannCurvatureTensor,
    SchoutenTensor,
    StressEnergyMomentumTensor,
    WeylTensor,
)

SYMS = symbols("t x y z")
T, X, Y, Z = SYMS
ETA = ImmutableMatrix(sympy.diag(-1, 1, 1, 1))


def boost(gamma, gamma_beta):
    # x' = L x, a boost along x
    return ImmutableMatrix(
        [
            [gamma, -gamma_beta, 0, 0],
            [-gamma_beta, gamma, 0, 0],
            [0, 0, 1, 0],
            [0, 0, 0, 1],
        ]
    )


# beta = 3/5, so gamma = 5/4 and gamma * beta = 3/4
L = boost(Rational(5, 4), Rational(3, 4))
L_INV = L.inv()
# old coordinates written in terms of the new ones
OLD_IN_NEW = dict(zip(SYMS, list(L_INV * sympy.Matrix(SYMS))))

# Curved metric whose t-x block is not boost invariant. It depends on y only,
# which the boost leaves alone, so the recomputation stays fast.
G_CURVED = ImmutableMatrix(sympy.diag(-(1 + Y**2), 1, 1, 1 + Y**2))
# Curved metric that depends on x, so the boost also mixes the coordinates
# inside the components.
G_MIXED = ImmutableMatrix(sympy.diag(-(1 + X**2), 1, 1 + X**2, 1))


def to_new_coords(arr):
    return sympy.Array(arr).subs(OLD_IN_NEW, simultaneous=True)


def assert_arrays_equal(a, b):
    diff = simplify(sympy.Array(a) - sympy.Array(b))
    assert all(e == 0 for e in sympy.flatten(diff))


@lru_cache(maxsize=None)
def metric_old(g):
    return MetricTensor(g.tolist(), SYMS, config="ll")


@lru_cache(maxsize=None)
def metric_new(g):
    # g'_ab = (L^-1)^c_a (L^-1)^d_b g_cd, written in the new coordinates
    g_new = simplify(to_new_coords(L_INV.T * g * L_INV))
    return MetricTensor(g_new.tolist(), SYMS, config="ll")


@lru_cache(maxsize=None)
def from_old_metric(cls, g=G_CURVED):
    return cls.from_metric(metric_old(g))


CURVATURE_CLASSES = [
    ChristoffelSymbols,
    RiemannCurvatureTensor,
    RicciTensor,
    EinsteinTensor,
    SchoutenTensor,
    WeylTensor,
]


def test_contraction_is_invariant():
    phi = symbols("phi")
    lt = boost(cosh(phi), sinh(phi)).tolist()
    u = GenericVector([1, 2, 3, 4], SYMS, config="u").lorentz_transform(lt)
    w = GenericVector([5, 6, 7, 8], SYMS, config="l").lorentz_transform(lt)
    assert simplify(sum(u[i] * w[i] for i in range(4))) == 70


def test_lowering_commutes_with_transform():
    eta = MetricTensor(ETA.tolist(), SYMS, config="ll")
    v = GenericVector([1, 2, 0, 0], SYMS, config="u", parent_metric=eta)
    expected = [Rational(1, 4), Rational(7, 4), 0, 0]
    lowered_then_boosted = v.change_config("l").lorentz_transform(L.tolist())
    boosted_then_lowered = v.lorentz_transform(L.tolist()).change_config("l")
    assert_arrays_equal(lowered_then_boosted.tensor(), expected)
    assert_arrays_equal(boosted_then_lowered.tensor(), expected)


@pytest.mark.parametrize("config", ["ll", "uu"])
def test_metric_matches_matrix_formula(config):
    g = metric_old(G_MIXED)
    if config == "uu":
        g = g.inv()
    boosted = g.lorentz_transform(L.tolist())
    assert boosted.config == config
    m = ImmutableMatrix(g.tensor().tolist())
    expected = L_INV.T * m * L_INV if config == "ll" else L * m * L.T
    assert_arrays_equal(boosted.tensor(), expected.tolist())


def test_boosted_perfect_fluid():
    # T_ab = (rho + p) u_a u_b + p eta_ab for a fluid at rest, seen from a
    # frame where it moves with u'^a = L^a_b u^b
    rho, p = symbols("rho p")
    eta = MetricTensor(ETA.tolist(), SYMS, config="ll")
    at_rest = sympy.diag(rho, p, p, p).tolist()
    fluid = StressEnergyMomentumTensor(at_rest, SYMS, config="ll", parent_metric=eta)
    boosted = fluid.lorentz_transform(L.tolist())
    u_new = ETA * L * ImmutableMatrix([1, 0, 0, 0])
    expected = (rho + p) * u_new * u_new.T + p * ETA
    assert isinstance(boosted, StressEnergyMomentumTensor)
    assert_arrays_equal(boosted.tensor(), expected.tolist())


@pytest.mark.parametrize("cls", CURVATURE_CLASSES)
def test_matches_quantity_from_boosted_metric(cls):
    boosted = from_old_metric(cls).lorentz_transform(L.tolist())
    assert isinstance(boosted, cls)
    recomputed = cls.from_metric(metric_new(G_CURVED))
    if recomputed.config != boosted.config:
        recomputed = recomputed.change_config(boosted.config)
    assert_arrays_equal(to_new_coords(boosted.tensor()), recomputed.tensor())


def test_christoffels_transform_as_tensor_when_coordinates_mix():
    # The inhomogeneous term in the Christoffel transformation law has
    # d^2 x' / dx dx, which is zero for a constant L
    boosted = from_old_metric(ChristoffelSymbols, G_MIXED).lorentz_transform(L.tolist())
    recomputed = ChristoffelSymbols.from_metric(metric_new(G_MIXED))
    assert_arrays_equal(to_new_coords(boosted.tensor()), recomputed.tensor())


@pytest.mark.parametrize(
    "cls", CURVATURE_CLASSES + [StressEnergyMomentumTensor, GenericVector]
)
def test_parent_metric_is_boosted(cls):
    if cls is GenericVector:
        obj = GenericVector([1, 0, 0, 0], SYMS, "u", parent_metric=metric_old(G_CURVED))
    else:
        obj = from_old_metric(cls)
    boosted = obj.lorentz_transform(L.tolist())
    assert boosted.parent_metric is not None
    assert_arrays_equal(
        boosted.parent_metric.tensor(),
        metric_old(G_CURVED).lorentz_transform(L.tolist()).tensor(),
    )


def test_change_config_commutes_with_transform():
    ch = from_old_metric(ChristoffelSymbols)
    a = ch.change_config("lll").lorentz_transform(L.tolist())
    b = ch.lorentz_transform(L.tolist()).change_config("lll")
    assert_arrays_equal(a.tensor(), b.tensor())


def test_ricci_scalar_is_invariant():
    ricci = from_old_metric(RicciTensor, G_MIXED).lorentz_transform(L.tolist())
    scalar = RicciScalar.from_riccitensor(ricci)
    original = RicciScalar.from_metric(metric_old(G_MIXED))
    assert simplify(scalar.expr - original.expr) == 0
