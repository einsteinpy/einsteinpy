from sympy import cos, sin, symbols

from einsteinpy.symbolic.metric import MetricTensor


def KerrNewmanDeSitter():
    """
    Kerr-Newman-(A)dS metric in natural units (G = M = c = k_e = 1).
    Coordinates: (t, r, theta, phi) in Boyer-Lindquist form.
    Signature: (-, +, +, +).

    Reference:
        Garnier, A. (2023).
        "Motion equations in a Kerr-Newman-de Sitter spacetime:
        some methods of integration and application to black holes shadowing in Scilab."
        Classical and Quantum Gravity, 40(13), 135011.
        doi:10.1088/1361-6382/accbfe

    """
    # Coordinates and parameters
    t, r, th, ph = symbols("t r theta phi", real=True)
    a, Q, Lambda = symbols("a Q Lambda", real=True)

    # Abbreviations
    lam = Lambda / 3
    chi = 1 + lam * a**2
    Sigma = r**2 + a**2 * cos(th) ** 2
    Delta_r = (1 - lam * r**2) * (r**2 + a**2) - 2 * r + Q**2  # M=1 in natural units
    Delta_th = 1 + lam * a**2 * cos(th) ** 2

    # Metric components g_{μν} in order (t, r, theta, phi)
    g_tt = (a**2 * sin(th) ** 2 * Delta_th - Delta_r) / (chi**2 * Sigma)
    g_rr = Sigma / Delta_r
    g_tht = Sigma / Delta_th
    g_tph = a * sin(th) ** 2 * (Delta_r - (r**2 + a**2) * Delta_th) / (chi**2 * Sigma)
    g_phph = (
        sin(th) ** 2
        * ((r**2 + a**2) ** 2 * Delta_th - a**2 * sin(th) ** 2 * Delta_r)
        / (chi**2 * Sigma)
    )

    metric = [
        [g_tt, 0, 0, g_tph],
        [0, g_rr, 0, 0],
        [0, 0, g_tht, 0],
        [g_tph, 0, 0, g_phph],
    ]

    coords = (t, r, th, ph)
    return MetricTensor(metric, coords, "ll", name="KerrNewmanDeSitterMetric")
