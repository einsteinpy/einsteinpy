import astropy.units as u
import numpy as np
from scipy.optimize import newton

from einsteinpy.ijit import jit


class Shadow:
    """
    Class for calculating shadow of black holes as seen by distant observer.
    """

    @u.quantity_input(
        mass=u.kg, distance=u.km,
    )
    def __init__(self, mass, n_rays, distance):
        self.mass = mass
        self.n_rays = n_rays
        self.distance = distance
        self.b_crit = 3 * np.sqrt(3) * self.mass
        self.b = self._compute_B()
        self.z = list()
        for i in self.b:
            root_eqn = lambda r: r / ((1 - (2 * self.mass.value / r))) ** 0.5 - i
            root = newton(root_eqn, 0.1)
            if np.isreal(root):
                self.z.append([i, np.real(root)])
        self.z = np.array(self.z)
        self.intensity = self._intensity()

    @jit
    def _compute_B(self):
        return np.linspace(self.b_crit.value, self.distance, self.n_rays)

    def _intensity_blue_sch(self, r, b):
        GTT = 1 - (2 * self.mass.value / r)
        GRR = (1 - (2 * self.mass.value / r)) ** (-1)
        KRKText = ((GTT / GRR) * (1 - (b ** 2 * GTT / (r ** 2)))) ** 0.5
        Gblue = (
            (1 / GTT) - KRKText * (GRR / GTT) * ((1 - GTT) / (GTT * GRR)) ** 0.5
        ) ** (-1)
        Iblue = -(Gblue ** 3) * (GTT / (r ** 2)) * (1 / KRKText)
        return Iblue

    def _intensity_red_sch(self, r, b):
        GTT = 1 - (2 * self.mass.value / r)
        GRR = (1 - (2 * self.mass.value / r)) ** (-1)
        KRKText = ((GTT / GRR) * (1 - (b ** 2 * GTT / (r ** 2)))) ** 0.5
        Gred = (
            (1 / GTT) + KRKText * (GRR / GTT) * ((1 - GTT) / (GTT * GRR)) ** 0.5
        ) ** (-1)
        Ired = (Gred ** 3) * (GTT / (r ** 2)) * (1 / KRKText)
        return Ired

    def _intensity(self):
        intensity = []
        for i in range(len(self.z)):
            b = self.z[i, 0]
            val1, err1 = si.integrate.quadrature(
                self._intensity_blue_sch, self.distance.value, self.z[i, 1], args=(b,)
            )
            val2, err2 = si.integrate.quadrature(
                self._intensity_red_sch, self.z[i, 1], self.distance.value, args=(b,)
            )
            intensity.append(val1 + val2)
        return intensity
