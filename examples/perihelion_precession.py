"""Perihelion precession in Schwarzschild spacetime.

Computes the geodesic of a test particle orbiting a central mass,
extracts perihelion positions from the trajectory, and measures the
apsidal precession angle.  The result is compared to the leading-order
Einstein formula.

Everything uses geometrized units (c = G = M = 1), so the
Schwarzschild radius is r_s = 2.

References
----------
- Einstein, A. (1915). "Explanation of the Perihelion Motion of Mercury
  from the General Theory of Relativity."
- MTW, *Gravitation* (1973), Section 25.6.

"""

import numpy as np

from einsteinpy.geodesic import Timelike


def find_perihelion_angles(phi, r, min_sep=5.0):
    """Return azimuthal angles at each perihelion (local minimum of *r*).

    Uses parabolic interpolation for sub-step precision and requires
    a minimum angular separation between successive detections to
    avoid double-dip artefacts.

    Parameters
    ----------
    phi : numpy.ndarray
        Azimuthal coordinate along the trajectory.
    r : numpy.ndarray
        Radial coordinate along the trajectory.
    min_sep : float, optional
        Minimum angular separation (radians) between successive
        perihelion detections. Default is 5.0 (just under 2*pi).

    Returns
    -------
    numpy.ndarray
        Angles (in radians) at which *r* is locally minimal.

    """
    angles = []
    for i in range(1, len(r) - 1):
        if r[i] < r[i - 1] and r[i] < r[i + 1]:
            # Parabolic interpolation for sub-step accuracy
            rm, r0, rp = r[i - 1], r[i], r[i + 1]
            denom = 2.0 * (rm - 2.0 * r0 + rp)
            if abs(denom) > 1e-12:
                frac = (rm - rp) / denom
                angle = phi[i] + frac * (phi[i + 1] - phi[i - 1]) / 2.0
            else:
                angle = phi[i]

            if len(angles) == 0 or (angle - angles[-1]) > min_sep:
                angles.append(angle)

    return np.array(angles)


def main():
    """Compute and verify perihelion precession in Schwarzschild spacetime."""
    # Moderately relativistic orbit: perihelion at r ~ 20 (= 10 r_s).
    # Angular momentum L = 5.5 gives eccentricity ~ 0.34.
    position = [20.0, np.pi / 2, 0.0]
    momentum = [0.0, 0.0, 5.5]

    geod = Timelike(
        metric="Schwarzschild",
        metric_params=(),
        position=position,
        momentum=momentum,
        steps=20000,
        delta=1.0,
        return_cartesian=False,
        suppress_warnings=True,
    )

    traj = geod.trajectory
    data = traj[1]

    r = data[:, 1]
    phi = data[:, 3]

    # Unwrap phi for monotonic growth
    phi = np.unwrap(phi)

    # Observe orbit apsides
    r_peri = r.min()
    r_apo = r.max()
    p = 2.0 * r_peri * r_apo / (r_peri + r_apo)
    e = (r_apo - r_peri) / (r_apo + r_peri)

    peri_angles = find_perihelion_angles(phi, r)

    if len(peri_angles) < 3:
        print("Not enough perihelion passages detected.")
        return

    # Drop first passage (transient) and use the rest
    peri_angles = peri_angles[1:]
    diffs = np.diff(peri_angles)
    precession_per_orbit = np.mean(diffs) - 2.0 * np.pi

    # Leading-order Einstein formula (geometrized: delta_phi = 6*pi*M / p).
    # This is the 1PN approximation; higher-order corrections of order
    # (M/p) ~ 4% enlarge the exact result by roughly 20% at this p.
    analytic_1pn = 6.0 * np.pi / p

    print("Perihelion precession in Schwarzschild spacetime")
    print("=" * 50)
    print("Geometrized units (c = G = M = 1, r_s = 2)")
    print()
    print("Orbit parameters (measured from trajectory):")
    print(f"  Perihelion radius  r_peri = {r_peri:.2f}")
    print(f"  Aphelion radius    r_apo  = {r_apo:.2f}")
    print(f"  Semi-latus rectum  p      = {p:.2f}")
    print(f"  Eccentricity       e      = {e:.4f}")
    print(f"  Angular momentum   L      = {momentum[2]}")
    print()
    print(f"Precession ({len(peri_angles)} perihelion passages, "
          f"{len(diffs)} intervals):")
    print(f"  Measured (numerical) : {precession_per_orbit:.6f} rad/orbit")
    print(f"  Einstein 1PN formula : {analytic_1pn:.6f} rad/orbit")
    rel_err = abs(precession_per_orbit - analytic_1pn) / analytic_1pn * 100
    print(f"  Relative difference  : {rel_err:.1f}%")
    print()
    print("Note: The ~20% difference is expected for an orbit this close")
    print("to the central mass (r_peri = 10 r_s). The 1PN formula is")
    print("only the leading-order term; higher-order post-Newtonian")
    print("corrections are significant at moderate field strengths.")
    print("The key result is that GR precession is clearly present")
    print("and prograde — Newtonian gravity predicts zero precession.")


if __name__ == "__main__":
    main()
