# Constant hyperresistive transport model

This page owns the equations and Nyquist heat convention for CAP-011 S2.
[The source-range plan](TRANSPORT_SOURCE_RANGE_PLAN.md) owns current implementation
status, canonical owners, compatibility, rejection and all acceptance gates.
Its composed implementation awaits GitHub validation. The public setup and supported
combinations are in the [transport guide](EVOLVING_ELECTRONS.md#constant-hyperresistive-transport).
The fixed-step first-order [ion-collision wrapper](COLLISIONAL_HYBRID.md) reuses
this model without an additional heat source or reservoir; its S2 admission and
essential regressions await exact-revision GitHub verification.

## Continuous and discrete balance

With J=∇×B, R=ηJ−ηH∇²J and periodic boundaries, magnetic work is
`−∫J·R = −η∫|J|² −ηH∫|∇J|²`. The continuum positive heat is
`H=η|J|²+ηH Σdc(∂d Jc)²`. The identity
`J·∇²J = ∇·Σc(Jc∇Jc)−Σdc(∂d Jc)²` implies the local total-energy
flux is the usual mechanical/enthalpy plus Poynting flux **minus**
`ηH Σc Jc∇Jc`. This is a dissipative energy-flux contribution; it does not add
an independent electron conduction closure to S2. Boundary integrals vanish
only for the supported periodic geometry.

Let Dd be the existing real skew-adjoint first derivative, and Pd the orthogonal
projection onto axis-d Nyquist modes (zero for odd axis lengths). Full L obeys
`−L=Σd(transpose(Dd)*Dd+kd,N² Pd)`, since Dd has multiplier ik away from Nyquist,
zero at Nyquist, and Pd²=Pd. Define Nd=kd,N Pd. Commuting real self-adjoint
projections preserve mixed modes independently along each axis. Therefore
`−<J,LJ>=Σdc(||Dd Jc||²+||Nd Jc||²)` with the grid-volume inner product.
S2 chooses cell heat `Hh=ηH Σdc[(Dd Jc)²+(Nd Jc)²]`. It is positive and
retains resolved-mode gradient heating; Nyquist heat is a stated collocation
convention, not a resolved continuum derivative. Real collocation product
rules and local Poynting identities do not hold exactly for arbitrary aliased
products; the accepted semidiscrete claim is the global SBP energy identity.
Time discretization remains RK4 with frozen moments and second-order complete
particle evolution, not exact fully discrete energy conservation.

For a mixed Nyquist magnetic mode, distinguish the real-curl wave vector κ
(zero on each Nyquist axis) from the full-Laplacian vector k. A transverse
mode has damping `ν=η|κ|²+ηH|κ|²|k|²`, bounded by
`η kmax²+ηH kmax⁴`. A pure curl-null Nyquist mode has J=0 and no dissipative
heat. The timestep neighbor gate checks this actual mixed-mode rate and the
production field/pressure RK4 polynomial against the parent's separate
combined-stability repair; it does not substitute `ηH|k|⁴` for the mixed rate.

Primary references checked before implementation:

- Stanier & Chacón, JCP 459 (2022), 111144,
  [equations (1)–(5)](https://arxiv.org/pdf/2110.03886): frictionless ion
  field, dissipative Faraday field and electron pressure/energy reservoirs.
  The earlier Stanier, Chacón & Chen JCP 376 (2019), 597–616
  [DOI metadata](https://doi.org/10.1016/j.jcp.2018.09.038) was checked;
  its full body was not independently retrieved for this slice.
- Stanier et al., Phys. Plasmas 20 (2013), 122302,
  [equations (3) and (5), section II.A](https://scientific-publications.ukaea.uk/wp-content/uploads/Published/POPvol20p122302.pdf):
  −ηH∇²J paired with positive ηH(∇J)² heating. This is a two-fluid reference;
  the Fourier/Nyquist identity above is independently derived here, not
  attributed to its finite-element discretization.
