import numpy as np

from einsteinpy.coordinates.utils import (
    bl_to_cartesian_fast,
    cartesian_to_bl_fast,
    cartesian_to_spherical_fast,
    spherical_to_cartesian_fast,
)


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
        return cartesian_to_spherical_fast(
            self.x_si,
            self.y_si,
            self.z_si,
            self.vx_si,
            self.vy_si,
            self.vz_si,
            self._velocities_provided,
        )

    def convert_bl(self, a):
        return cartesian_to_bl_fast(
            self.x_si,
            self.y_si,
            self.z_si,
            a,
            self.vx_si,
            self.vy_si,
            self.vz_si,
            self._velocities_provided,
        )


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
        return spherical_to_cartesian_fast(
            self.r_si,
            self.t_si,
            self.p_si,
            self.vr_si,
            self.vt_si,
            self.vp_si,
            self._velocities_provided,
        )

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
        return bl_to_cartesian_fast(
            self.r_si,
            self.t_si,
            self.p_si,
            self.a_si,
            self.vr_si,
            self.vt_si,
            self.vp_si,
            self._velocities_provided,
        )

    def convert_spherical(self):
        transformed_cartesian = self.convert_cartesian()
        cart = CartesianConversion(*transformed_cartesian)
        return cart.convert_spherical()
