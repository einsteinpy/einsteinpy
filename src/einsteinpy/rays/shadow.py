import warnings

import numpy as np
from astropy import units as u
from scipy.integrate import fixed_quad
from scipy.interpolate import interp1d
from scipy.optimize import newton


class Shadow:
    """
    Class for plotting the shadow of Schwarzschild Black Hole surrounded by a
    thin accreting emission disk as seen by a distant observer.
    """

    @u.quantity_input(mass=u.kg, fov=u.km)
    def __init__(self, mass, n_rays, fov, limit=0.001):
        self.mass = mass.to(u.kg)
        self.limit = limit
        self.n_rays = n_rays
        self.fov = fov.to(u.km)
        self.horizon = 2 * self.mass.value  # To be changed after 0.3.0
        self.b_crit = 3 * np.sqrt(3) * self.mass
        self.b = self._compute_B()
        self.z = list()
        self.bfin = list()
        for i in self.b:
            root = newton(self._root_equation, 0.1, args=(i,))
            if np.isreal(root):
                self.bfin.append(i)
                self.z.append([i, np.real(root)])
        self.z = np.array(self.z)
        self.k0 = self._intensity()
        self.k1 = self._intensity_from_event_horizon()
        self.intensity = self.k1 + self.k0
        # Just to make the plot symmetric on -x axis
        self.fb1 = list(self.b2) + list(self.bfin)
        self.fb2 = np.asarray(list(-np.asarray(self.b2)) + list(-np.asarray(self.bfin)))

    def _compute_B(self):
        """
        Returns an array of impact parameters
        """
        return np.linspace(self.b_crit.value, self.fov.value, self.n_rays)

    def _root_equation(self, r_tp, i):
        """
        Returns the root of the equation for ``r_tp`` (turning points) for some impact parameter
        """
        # emath.sqrt is domain-agnostic and this is a complex equation
        return r_tp / np.emath.sqrt(1 - (2 * int(self.mass.value) / r_tp)) - i

    def _intensity_blue_sch(self, r, b):
        """
        Returns the integrand for the blue shifted intensity to be integrated.
        Reference : Cosimo Bambi, 10.1103/PhysRevD.87.107501
        """
        GTT = 1 - (2 * self.mass.value / r)
        GRR = (1 - (2 * self.mass.value / r)) ** (-1)
        KRKText = ((GTT / GRR) * (1 - (b**2 * GTT / (r**2)))) ** 0.5
        Gblue = (
            (1 / GTT) - KRKText * (GRR / GTT) * ((1 - GTT) / (GTT * GRR)) ** 0.5
        ) ** (-1)
        Iblue = -(Gblue**3) * (GTT / (r**2)) * (1 / KRKText)
        return Iblue

    def _intensity_red_sch(self, r, b):
        """
        Returns the integrand for the red shifted intensity to be integrated.
        Reference : Cosimo Bambi, 10.1103/PhysRevD.87.107501
        """
        GTT = 1 - (2 * self.mass.value / r)
        GRR = (1 - (2 * self.mass.value / r)) ** (-1)
        KRKText = ((GTT / GRR) * (1 - (b**2 * GTT / (r**2)))) ** 0.5
        Gred = (
            (1 / GTT) + KRKText * (GRR / GTT) * ((1 - GTT) / (GTT * GRR)) ** 0.5
        ) ** (-1)
        Ired = (Gred**3) * (GTT / (r**2)) * (1 / KRKText)
        return Ired

    def _intensity(self):
        """
        Returns an array of the integrated values using ~scipy.integrate.quadrature as the
        intensities for the blue shifted and red shifted rays above the critical impact paratmeter
        from the distance to the emitter
        """
        intensity = []
        for i in np.arange(len(self.z)):
            b = self.z[i][0]
            val1, _ = fixed_quad(
                self._intensity_blue_sch, self.fov.value, self.z[i][1], args=(b,)
            )
            val2, _ = fixed_quad(
                self._intensity_red_sch, self.z[i][1], self.fov.value, args=(b,)
            )
            intensity.append(val1 + val2)
        return intensity

    def _intensity_from_event_horizon(self):
        """
        Returns an array of the integrated values using ~scipy.integrate.quadrature as the
        intensities for the blue shifted and red shifted rays below the critical impact paratmeter
        from the event horizon to the distance given.
        """
        self.b2 = np.linspace(self.limit, self.b_crit.value, len(self.bfin))
        k1 = list()
        for i in self.b2:
            arg = i
            val3, _ = fixed_quad(
                self._intensity_red_sch, self.horizon, self.fov.value, args=(arg,)
            )
            k1.append(val3)
        return k1

    def smoothen(self, points=500):
        """
        Sets the interpolated values for the intensities for smoothening of the plot
        using ~scipy.interpolate.interp1d
        """
        b_new = np.linspace(np.min(self.fb1), np.max(self.fb1), points)
        interpolation = interp1d(self.fb1, self.intensity, kind="cubic")
        smoothened = interpolation(b_new)
        self.intensity = smoothened
        self.fb1 = b_new
        self.fb2 = -1 * b_new
