# Kerr–Newman charge scaling and Lorentz force

This note records the unit conventions, the source of the charge-scaling problem in issue [#144](https://github.com/einsteinpy/einsteinpy/issues/144), the implementation changes in [PR #687](https://github.com/einsteinpy/einsteinpy/pull/687) and its follow-up revision, and the remaining caveats.

## Summary

The numerical factor reported in issue #144 is connected to converting a specific charge between SI and geometrized conventions. It is not a dimensionless correction factor that can be applied indiscriminately to either charge parameter.

The implementation changes address several distinct problems:

1. The Lorentz-force term must contract one index of the field tensor with the inverse metric before acting on the contravariant four-velocity.
2. The Kerr–Newman metric must include the source charge through the \(r_Q^2\) term in \(\Delta\).
3. The listed electromagnetic tensor components must agree with \(F=dA\) for the implemented potential.
4. The field used in the SI-coordinate equation of motion must be converted consistently from the geometrized field convention.

The electromagnetic potential and field must preserve the sign of \(Q\), while the metric depends only on \(Q^2\). The follow-up implementation uses signed \(r_Q=Q\sqrt{Gk_e}/c^2\) in the potential and field and retains \(r_Q^2\) in the metric. Reversal tests and opposite-sign force tests protect that distinction.

## Conventions used by the metric class

The metric class accepts SI quantities: mass \(M\) in kilograms, source charge \(Q\) in coulombs, and Boyer–Lindquist radius \(r\) in metres. Its coordinate vector uses \((t,r,\theta,\phi)\), with \(t\) in seconds. This is why the metric components include \(g_{tt}\propto c^2\) and \(g_{t\phi}\propto c\).

The standard geometrized Kerr–Newman expressions are commonly written with \(x^0=ct\), and with \(G=c=k_e=1\), where \(k_e=1/(4\pi\epsilon_0)\). In those expressions, mass and charge are represented by lengths. For SI inputs, the relevant scales are

\[
r_s=\frac{2GM}{c^2},\qquad
\alpha=\frac{GMa}{c^2},\qquad
r_Q=\frac{Q\sqrt{Gk_e}}{c^2}.
\]

The metric uses

\[
\Delta=r^2-r_s r+\alpha^2+r_Q^2.
\]

Only \(r_Q^2\) belongs in the metric. The electromagnetic potential and field are linear in the source charge and therefore require the sign of \(Q\) to be preserved.

## Specific-charge conversion and the reported factor

For a source mass \(M\), the geometrized mass length is \(GM/c^2\), while the geometrized charge length is \(Q\sqrt{Gk_e}/c^2\). Their dimensionless ratio is

\[
\frac{r_Q}{GM/c^2}=\frac{Q}{M}\sqrt{\frac{k_e}{G}}.
\]

The same factor converts a test particle's SI specific charge to its dimensionless geometrized value:

\[
\left(\frac{q}{m}\right)_{\rm geom}
=\left(\frac{q}{m}\right)_{\rm SI}\sqrt{\frac{k_e}{G}}.
\]

Here \(\sqrt{k_e/G}\) has units of kg/C; its numerical value is about \(1.16\times10^{10}\) for the constants used when the original value in issue #144 was recorded. It is not dimensionless. In particular, converting a charge \(Q\) in coulombs to a geometric length requires \(\sqrt{Gk_e}/c^2\), not multiplication by \(\sqrt{k_e/G}\) alone. The factor in issue #144 is evidence of a mixed unit convention in the force calculation, not a universal charge fudge factor.

## Lorentz force in the SI-coordinate equation

For a particle with SI specific charge \(\eta=(q/m)_{\rm SI}\), the contravariant equation of motion is

\[
\frac{d^2x^\mu}{d\tau^2}
=-\Gamma^\mu{}_{\alpha\beta}u^\alpha u^\beta
+\eta\,g^{\mu\alpha}F_{\alpha\nu}u^\nu,
\qquad u^\mu=\frac{dx^\mu}{d\tau}.
\]

The contraction \(g^{\mu\alpha}F_{\alpha\nu}u^\nu\) is essential: it raises the first index of the covariant field tensor and leaves the second index available to contract with the contravariant four-velocity. A bare matrix product \(F^{\mu\nu}u^\nu\) repeats an upper index and is not the required tensor contraction. Appending the covariant metric to that product does not repair the index structure.

The electromagnetic expressions are naturally written in the geometrized \(x^0=ct\) convention, while the metric class's ODE state uses \(t\) in seconds. Pulling a covector or two-form back to the \(t\) coordinate contributes a factor of \(c\) for each time index. In addition, converting the geometrized field scale for use with SI \(q/m\) contributes the specific-charge conversion above and the factor of \(c\) associated with the proper-time parametrization. The follow-up implementation applies these factors before evaluating the mixed-index Lorentz term.

This conversion should be kept together with its derivation in code comments and tests. The conversion of a time-index component and the overall unit conversion are separate operations; they should not be collapsed into an unexplained numeric constant.

## Electromagnetic field consistency

The Kerr–Newman potential implemented by the class has the geometrized form

\[
A=\frac{r r_Q}{\rho^2}\left(dx^0-a\sin^2\theta\,d\phi\right),
\qquad \rho^2=r^2+a^2\cos^2\theta,
\]

up to the code's chosen overall sign convention. The covariant field tensor must be its exterior derivative, \(F=dA\), using that same convention. The follow-up revision corrects the signs of \(F_{r\phi}\) and \(F_{\theta\phi}\) so they agree with the implemented potential. A derivative-based regression test checks this relation.

The potential and every field component reverse sign under \(Q\mapsto -Q\). The metric remains unchanged under that transformation. Regression tests check both properties and, in the weak-field Reissner–Nordström limit, verify the Lorentz-force direction for all source/test charge sign combinations.

## Changes and regression checks

The implementation represented by the follow-up revision:

- passes the source charge to \(\Delta\) in both the metric and its radial derivatives;
- forms the Lorentz term as \(g^{-1}F u\), after the field conversion;
- corrects the two field-component signs that disagreed with \(F=dA\);
- adds checks for the Reissner–Nordström metric, the exterior derivative of the potential, static charged-particle equilibrium, and conservation of the stationary and axial canonical momenta along an integrated orbit.

Those tests are useful, but most of their relativistic reference values are expressed using the same geometrization relations as the implementation. Add an independent weak-field SI check so that a shared conversion error cannot pass on both sides of the comparison. For example, at rest at \(r=10^3\,\mathrm m\), use \(Q=+1\,\mathrm C\), \(q/m=+1\,\mathrm{C/kg}\), and a source mass small enough that gravity is negligible. The radial coordinate acceleration should approach

\[
\ddot r\simeq \frac{k_e Q(q/m)}{r^2}
\]

or about \(+8.99\times10^3\,\mathrm{m/s^2}\). Repeating with \(Q=-1\,\mathrm C\) should reverse the sign. A relative tolerance near \(10^{-6}\) is reasonable in this weak-field setup, but the test should use the library's own Coulomb constant and should confirm the initial state is normalized for proper-time evolution.

The canonical momentum components are conserved in a stationary, axisymmetric Kerr–Newman field. With the usual energy convention, \(E=-\pi_t\), while the axial angular momentum is \(L=\pi_\phi\). A test that checks constancy of \(\pi_t\) is valid; its comment should call it the conserved time component or use \(-\pi_t\) when naming it energy.

## Public electromagnetic-method unit contract

The force implementation converts the field locally in the ODE, but the public methods `em_potential_covariant` and `em_tensor_covariant` continue to return the geometrized expressions. Their docstrings should say which coordinate convention and units they use.

The contravariant methods currently obtain their outputs by raising those covariant results with the metric class's SI-coordinate inverse metric. That combines conventions and does not define a consistent geometrized or SI contravariant tensor. Documenting the covariant methods alone does not fix this. The contravariant methods should either convert to a consistent SI-coordinate covariant quantity before raising, raise using a metric in the geometrized convention, or be revised/deprecated so callers are not offered a hybrid result as a physical tensor.

The Kerr–Newman call sites normalize `Q` with `to_value(u.C)` before passing it to the metric helper or using it in \(r_Q\). This matters because `BaseMetric.delta` accepts `Quantity` objects through `primitive`, but `primitive` strips `.value` without converting to coulombs; a direct call with `1 * u.mC` would otherwise be treated as `1 C`. The same broader review should be applied to mass inputs, since `.value` is used elsewhere with an implicit kilogram assumption.

## References

- Issue [#144](https://github.com/einsteinpy/einsteinpy/issues/144): reported charge-scaling problem.
- [PR #687](https://github.com/einsteinpy/einsteinpy/pull/687) and its follow-up revision: implementation changes.
- Hackmann, E. and Xu, H., [“Charged particle motion in Kerr–Newman space-times”](https://arxiv.org/abs/1304.2142), *Physical Review D* 87, 124030 (2013).
- Liang, J. and Zhou, S., charged-particle equations and conserved quantities in Kerr–Newman spacetime, summarized in [“On neutralization of charged black holes”](https://academic.oup.com/mnras/article/488/2/2722/5531777).
- The extremal Reissner–Nordström / Majumdar–Papapetrou equilibrium is discussed in [“Charged massive particle at rest in the field of a Reissner–Nordström black hole”](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.75.044012).
