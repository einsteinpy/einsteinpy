from sympy import *
from einsteinpy.symbolic.tensor import *
from einsteinpy.symbolic.metric import *
from einsteinpy.symbolic.partial import *
init_printing()

t,r,th,ph = symbols('t r theta phi')
coords = [t, r, th, ph]
schw = diag(1-1/r, -1/(1-1/r), -r**2, -r**2*sin(th)**2)
g = Metric('g', coords, schw)
x = Tensor('x', coords, g)
y = Tensor('y', tensorproduct(coords, coords), g)
mu, nu, sigma, rho, lamda = indices('mu nu sigma rho lambda', g)
d = g.partial
gamma = g.christoffel
R = g.riemann
# RR = g.ricci_tensor
# RRR = g.ricci_scalar
# C = g.weyl
