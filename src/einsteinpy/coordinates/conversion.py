import numpy as np


class CartesianConversion:
    def __init__(self, x, y, z, v_x=None, v_y=None, v_z=None):
        self.x_si = x
        self.y_si = y
        self.z_si = z
        self.vx_si = v_x
        self.vy_si = v_y
        self.vz_si = v_z
        self._velocities_provided = not (
            (v_x is None) or (v_y is None) or (v_z is None)
        )

    def values(self):
        if self._velocities_provided:
            return self.x_si, self.y_si, self.z_si, self.vx_si, self.vy_si, self.vz_si
        return self.x_si, self.y_si, self.z_si

    def convert_spherical(self):
        r = np.sqrt(self.x_si ** 2 + self.y_si ** 2 + self.z_si ** 2)
        theta = np.arccos(self.z_si / r)
        phi = np.arctan2(self.y_si, self.x_si)
        if self._velocities_provided:
            n1 = self.x_si ** 2 + self.y_si ** 2
            n2 = n1 + self.z_si ** 2
            v_r = (
                self.x_si * self.vx_si + self.y_si * self.vy_si + self.z_si * self.vz_si
            ) / np.sqrt(n2)
            v_t = (
                self.z_si * (self.x_si * self.vx_si + self.y_si * self.vy_si)
                - n1 * self.vz_si
            ) / (n2 * np.sqrt(n1))
            v_p = -1 * (self.vx_si * self.y_si - self.x_si * self.vy_si) / n1

            return r, theta, phi, v_r, v_t, v_p
        return r, theta, phi

    def convert_bl(self, a):
        w = (self.x_si ** 2 + self.y_si ** 2 + self.z_si ** 2) - (a ** 2)
        r = np.sqrt(0.5 * (w + np.sqrt((w ** 2) + (4 * (a ** 2) * (self.z_si ** 2)))))
        theta = np.arccos(self.z_si / r)
        phi = np.arctan2(self.y_si, self.x_si)
        if self._velocities_provided:
            dw_dt = 2 * (
                self.x_si * self.vx_si + self.y_si * self.vy_si + self.z_si * self.vz_si
            )
            v_r = (1 / (2 * r)) * (
                (dw_dt / 2)
                + (
                    (w * dw_dt + 4 * (a ** 2) * self.z_si * self.vz_si)
                    / (2 * np.sqrt((w ** 2) + (4 * (a ** 2) * (self.z_si ** 2))))
                )
            )
            v_t = (-1 / np.sqrt(1 - np.square(self.z_si / r))) * (
                (self.vz_si * r - v_r * self.z_si) / (r ** 2)
            )
            v_p = (1 / (1 + np.square(self.y_si / self.x_si))) * (
                (self.vy_si * self.x_si - self.vx_si * self.y_si) / (self.x_si ** 2)
            )

            return r, theta, phi, v_r, v_t, v_p, a
        return r, theta, phi, a


class SphericalConversion:
    def __init__(self, r, theta, phi, v_r=None, v_t=None, v_p=None):
        self.r_si = r
        self.t_si = theta
        self.p_si = phi
        self.vr_si = v_r
        self.vt_si = v_t
        self.vp_si = v_p
        self._velocities_provided = not (
            (v_r is None) or (v_t is None) or (v_p is None)
        )

    def values(self):
        if self._velocities_provided:
            return self.r_si, self.t_si, self.p_si, self.vr_si, self.vt_si, self.vp_si
        return self.r_si, self.t_si, self.p_si

    def convert_cartesian(self):
        x = self.r_si * np.cos(self.p_si) * np.sin(self.t_si)
        y = self.r_si * np.sin(self.p_si) * np.sin(self.t_si)
        z = self.r_si * np.cos(self.t_si)
        if self._velocities_provided:
            v_x = (
                np.sin(self.t_si) * np.cos(self.p_si) * self.vr_si
                - self.r_si * np.sin(self.t_si) * np.sin(self.p_si) * self.vp_si
                + self.r_si * np.cos(self.t_si) * np.cos(self.p_si) * self.vt_si
            )
            v_y = (
                np.sin(self.t_si) * np.sin(self.p_si) * self.vr_si
                + self.r_si * np.cos(self.t_si) * np.sin(self.p_si) * self.vt_si
                + self.r_si * np.sin(self.t_si) * np.cos(self.p_si) * self.vp_si
            )
            v_z = (
                np.cos(self.t_si) * self.vr_si
                - self.r_si * np.sin(self.t_si) * self.vt_si
            )
            return x, y, z, v_x, v_y, v_z
        return x, y, z

    def convert_bl(self, a):
        transformed_cartesian = self.convert_cartesian()
        cart = CartesianConversion(*transformed_cartesian)
        return cart.convert_bl(a)


class BoyerLindquistConversion:
    def __init__(self, r, theta, phi, v_r=None, v_t=None, v_p=None, a=0.0):
        self.r_si = r
        self.t_si = theta
        self.p_si = phi
        self.vr_si = v_r
        self.vt_si = v_t
        self.vp_si = v_p
        self.a_si = a
        self._velocities_provided = not (
            (v_r is None) or (v_t is None) or (v_p is None)
        )

    def values(self):
        if self._velocities_provided:
            return (
                self.r_si,
                self.t_si,
                self.p_si,
                self.vr_si,
                self.vt_si,
                self.vp_si,
                self.a_si,
            )
        return self.r_si, self.t_si, self.p_si, self.a_si

    def convert_cartesian(self):
        xa = np.sqrt(self.r_si ** 2 + self.a_si ** 2)
        sin_norm = xa * np.sin(self.t_si)
        x = sin_norm * np.cos(self.p_si)
        y = sin_norm * np.sin(self.p_si)
        z = self.r_si * np.cos(self.t_si)
        if self._velocities_provided:
            v_x = (
                (self.r_si * self.vr_si * np.sin(self.t_si) * np.cos(self.p_si) / xa)
                + (xa * np.cos(self.t_si) * np.cos(self.p_si) * self.vt_si)
                - (xa * np.sin(self.t_si) * np.sin(self.p_si) * self.vp_si)
            )
            v_y = (
                (self.r_si * self.vr_si * np.sin(self.t_si) * np.sin(self.p_si) / xa)
                + (xa * np.cos(self.t_si) * np.sin(self.p_si) * self.vt_si)
                + (xa * np.sin(self.t_si) * np.cos(self.p_si) * self.vp_si)
            )
            v_z = (self.vr_si * np.cos(self.t_si)) - (
                self.r_si * np.sin(self.t_si) * self.vt_si
            )
            return x, y, z, v_x, v_y, v_z
        return x, y, z

    def convert_spherical(self):
        transformed_cartesian = self.convert_cartesian()
        cart = CartesianConversion(*transformed_cartesian)
        return cart.convert_spherical()
