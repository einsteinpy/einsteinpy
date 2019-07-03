from sympy import *

from einsteinpy.symbolic.metric import *
from einsteinpy.symbolic.partial import *
from einsteinpy.symbolic.tensor import *

init_printing()

t, r, th, ph, x, y, z = symbols("t r theta phi x y z", real=True)
momentum = symbols("E p_1:4", positive=True)
E, p1, p2, p3 = momentum
# schw = diag(1-1/r, -1/(1-1/r), -r**2, -r**2*sin(th)**2)
alpha, beta = symbols("alpha beta", cls=Function)
schw = diag(-exp(2 * alpha(r)), exp(2 * beta(r)), r ** 2, r ** 2 * sin(th) ** 2)
g = Metric("g", [t, r, th, ph], schw)
x = Tensor("x", [t, r, th, ph], g)
p = Tensor("p", momentum, g)
mu, nu, si, rh, la = indices("mu nu sigma rho lambda", g)
d = g.partial
Gamma = g.christoffel
R = g.riemann
G = g.einstein
C = g.weyl

# mink = diag(1, -1, -1, -1)
# eta = Metric('eta', [t, x, y, z], mink)
# al, be, ga = indices('alpha beta gamma', eta)
# E1, E2, E3, B1, B2, B3 = symbols('E_1:4 B_1:4', real=True)
# matrix = [[0, -E1, -E2, -E3], [E1, 0, -B3, B2], [E2, B3, 0, -B1], [E3, -B2, B1, 0]]
# F = Tensor('F', matrix, eta, symmetry=[[2]])
