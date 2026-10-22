"""Base classes for coordinate representations.

These provide shared logic for position coordinates and differentials,
reducing duplication across Cartesian, Spherical, and Boyer-Lindquist variants.
See https://github.com/einsteinpy/einsteinpy/issues/650
"""

from typing import Any, Dict, Tuple

import numpy as np
from astropy import units as u

from einsteinpy import constant
from einsteinpy.coordinates.utils import v0
from einsteinpy.utils import CoordinateError

_c = constant.c.value


class BaseSpatialCoordinate:
    """
    Mixin providing common behavior for position-only coordinate classes.

    Expects subclasses to set ``_dimension``, ``_dimension_order``, and ``t``
    before use. The first element of ``_dimension_order`` is the time component;
    the rest are spatial components in order.
    """

    _dimension: Dict[str, Any]
    _dimension_order: Tuple[str, ...]
    t: Any  # astropy Quantity

    def __getitem__(self, item: Any) -> Any:
        """
        Return coordinate components by name or integer index.

        Parameters
        ----------
        item : str or int
            Parameter name (e.g. ``"x"``, ``"r"``) or index (0 = t, 1–3 = spatial).
            Negative indices are supported.

        Returns
        -------
        ~astropy.units.Quantity or str
            The requested component, or ``system`` when ``"system"`` is passed.
        """
        if isinstance(item, (int, np.integer)):
            return self._dimension[self._dimension_order[item]]
        return self._dimension[item]

    def position(self) -> Tuple[float, ...]:
        """
        Return the position 4-vector in SI units (c*t, spatial...).

        Returns
        -------
        tuple
            4-tuple containing the position 4-vector in SI units.
        """
        t_val = _c * self.t.si.value
        spatial_names = self._dimension_order[1:]
        spatial_vals = [self._dimension[n].si.value for n in spatial_names]
        return (t_val, *spatial_vals)


class BaseCoordinateDifferential:
    """
    Mixin providing common behavior for coordinate differentials (position + velocity).

    Subclasses must set ``_spatial_component_names`` (class attr) and implement
    ``_velocity_si_values()`` returning (v1, v2, v3) in SI units.
    """

    _spatial_component_names: Tuple[str, ...] = ()  # override in subclass
    _v_t: Any
    t: Any  # astropy Quantity
    system: str

    def _velocity_si_values(self) -> Tuple[float, float, float]:
        """Return (v1, v2, v3) in SI units. Subclasses must implement."""
        raise NotImplementedError

    def position(self) -> Tuple[float, ...]:
        """Return the position 4-vector in SI units (c*t, spatial...)."""
        t_val = _c * self.t.si.value
        spatial_vals = [
            getattr(self, n).si.value for n in self._spatial_component_names
        ]
        return (t_val, *spatial_vals)

    @property
    def v_t(self) -> Any:
        """Return the timelike component of 4-velocity."""
        return self._v_t

    @v_t.setter
    def v_t(self, args: Any) -> None:
        g = args[0]
        if self.system != g.coords.system:
            raise CoordinateError(
                f"Metric object has been instantiated with a coordinate system, "
                f"({g.coords.system}) other than {self.system} Coordinates."
            )
        g_cov_mat = g.metric_covariant(self.position())
        v_t_val = v0(g_cov_mat, *self._velocity_si_values())
        self._v_t = v_t_val * u.m / u.s

    def velocity(self, metric: Any) -> Tuple[float, ...]:
        """Return the velocity 4-vector in SI units."""
        self.v_t = (metric,)
        v1, v2, v3 = self._velocity_si_values()
        return (self._v_t.value, v1, v2, v3)
