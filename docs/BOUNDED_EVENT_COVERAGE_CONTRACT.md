# Carried-error contract for one event cluster

Reviewed 2026-11-02. **Conditional contract; implemented source awaits
exact-revision GitHub acceptance.**
[`YeeImpactLongitudinalPIC`](YEE_IMPACT_LONGITUDINAL.md) implements these carried
checks for zero contacts or one connected reflection cluster per external step
on a closed longitudinal 1D grid with finite clouds and nonrelativistic CPU
Float32/64 particles. The derivation and source share the canonical cloud/current, original-force, exact source-row and
work owners. Source review does not establish runtime correctness, convergence,
JET/compiler behavior or complete-step resource costs.

**Disposition:** under the acceptance premises below, this carried-error algebra certifies physical no-hidden-contact coverage from the original external-step stored state through its approximate cluster and full remaining duration. It never restarts the physical proof with zero errors at the numerical event. Conservative or inconclusive checks reject. This is a conditional mathematical contract, not proof that an implementation enforces it, runtime acceptance, or a global entire-run initial-trajectory/convergence theorem.

## 1. Reference and independent global bounds

The reference is the autonomous fixed-source finite-width longitudinal semidiscretization initialized from exact rational interpretations of this external step's stored x0/v0/E0. Set `Qp=q_p*w_p`, `alpha_p=q_p/m_p`, `h=Q(L)/N`, `delta=Q(stored_dx)`, and `a_cloud` to the original positive physical half-width. Plates/background/sources/counts remain fixed; reflection changes only velocity sign at the center's wall. There is no absorption, collision, antenna or explicit time-dependent forcing. Tau below is the resolved stored external duration interpreted exactly, rather than a reconstructed difference of public floating clocks.

From the [triangular cloud primitives](../src/coupling/shapes.jl) and the
[canonical endpoint-current and field-gather owners](../src/coupling/yee_charge.jl),

```
0<=g_e<=M=min(1/h,1/a_cloud)
sum g<=1/h,     sum |g'_e|<=2/(a_cloud*h)
R_e=E_e+sum Qp*A_e(xp)             (exact physical invariant).
```

For the entire external horizon tau use the already independently derived bounds

```
E*=norm(E0,Inf)+sum |Qp|
A_p=|alpha_p|*delta*E*/h
V_p=|v0p|+A_p*tau
J_p=|alpha_p|*delta*[(M/h)*sum_j |Qj|*V_j
                                   +2*E*/(a_cloud*h)*V_p].
```

They bound physical field, acceleration, speed and the Lipschitz constant of acceleration until a particle's own next hit, including across other particles' reflections. They are formed from the external state, not from the desired remainder conclusion. The field invariant bounds E. The cloud bounds above give
`E_dot,e=-sum_j Qj*g_e(xj)*vj` and
`a_dot,p=alpha_p*delta*sum_e(E_dot,e*g_e(xp)+E_e*g_prime,e(xp)*vp)`.
Their absolute sums give J. Reflections preserve speed magnitude; other
particles' hits can jump the force derivative while keeping force continuous.
These Lipschitz bounds do not require a rest-mode assumption.

## 2. First-hit/hull premises and within-cluster coverage

First independently certify a nonempty selected set C. Every selected particle has one assigned wall, its **first physical hit** `t_i in [l_i,u_i]`, and incoming physical speed at least `nu_hit,i>0`. The physical-bracket, first-prefix, nongrazing, horizon and branch checks must provide these actual certificates. Every nonselected particle must have no first hit on either wall through the hull upper time b. Define

```
a=min_i l_i,   b=max_i u_i,   W=b-a
0<a<=b<tau,    s in [a,b]
Delta=max(s-a,b-s),       bn=tau-s>0.
```

Here a is a time, distinct from a_cloud. The chosen earliest predictor upper endpoint lies in its own widened physical interval and therefore in the hull; this inclusion must be checked, not inferred from a failed/clipped bracket. Physical intervals establish earliest cluster ordering and exclude unresolved earlier/omitted contacts. Each particle has at most one selected side. Complete all numerical first-chunk kicks/coordinates/shared field/current at s; selected x1 equals its assigned wall exactly, nonselected x1 is strict interior, and each selected completed numerical vminus passes its incoming-sign gate. Put vplus=-vminus for C, vplus=vminus otherwise. The force a1 is unchanged by this velocity-only sign flip.

Additionally require, for every selected particle,

```
A_i*W < nu_hit,i/2,       V_i*W < L.
```

After its own true hit, inward speed stays above `nu_hit,i-A_i*W>0` through b, excluding return to that wall. Its traveled distance is below L, excluding the opposite wall. Apply these estimates until any hypothetical next own hit to contradict such a hit; no second-hit absence is assumed to prove itself. With nonselected exclusion through b, these gates establish the common true postcluster time `t_g=max_i t_i in [a,b]`, with one hit for each selected particle and none for others through t_g (indeed through b). The physical hits need not be simultaneous. No event-loop or Zeno claim enters this bounded proof.

## 3. Measured map/source residuals and uniform hull errors

Let x1/E1/vminus be the actual completed first numerical chunk at s, before numerical reflection; a0 and a1 are canonical exact original-factor gathers at external start and x1/E1. For all particles measure in Q

```
dx_p=Q(x1p)-[Q(x0p)+s*Q(v0p)+s^2*a0p/2]
kappa_p=Q(vminus_p)-Q(v0p)-s*(a0p+a1p)/2.
```

These include snap, coordinate and both kick rounding errors. Obtain the unrounded **canonical exact-kernel** endpoint impulse

```
Istar_e=sum_p Qp*[A_e(x1p)-A_e(x0p)]
rho_e=Q(E1_e)-Q(E0_e)+Istar_e.
```

Rho is the actual numerical R-invariant defect, including current and field-entry rounding. Q of rounded current or Float64 workspace sums is not Istar. Use actual s/endpoints, not a trial root's residuals.

Uniformly for true physical times t in [a,b], set

```
X_p=V_p*W                                     (selected)
X_p=|dx_p|+J_p*s^3/6+V_p*Delta                (nonselected).
```

A selected true position equals its numerical wall at its own t_i; bounded speed gives distance <=V W to that wall at any hull time, before or after reflection. A nonselected particle has no own hit through b: its initial Taylor tube at s gives `|dx|+J*s^3/6`, and motion to any hull time adds V Delta. Thus `|xphys,p(t)-x1p|<=X_p` without pretending the true hits share one time.

For a nonnegative vector X define

```
L_p(X)=|alpha_p|*delta*[(M/h)*sum_j |Qj|*X_j
                                  +2*E*/(a_cloud*h)*X_p]
B_R,p=|alpha_p|*delta*norm(rho,Inf)/h
D_a,p=B_R,p+L_p(X).
```

The exact invariant gives `norm(E1-Ephys(t),Inf)<=norm(rho,Inf)+M*sum |Qj|X_j`. In the force difference use

```
(E1-Ephys)*g(x1)+Ephys*[g(x1)-g(xphys)].
```

The second coefficient is the physical Ephys bounded by E*, not an unproved bound on E1. Basis sums then prove `|a1p-aphys,p(t)|<=D_a,p` uniformly on the hull; no cross term is missing.

## 4. Incoming velocity, reflection and common true postcluster state

For a Lipschitz acceleration before that particle's own first hit,

```
int_0^t a(u)du-t*[a(0)+a(t)]/2
       =int_0^t (t/2-u)*a'(u)du
absolute error <=J_p*t^2/4.
```

The derivative exists almost everywhere; other particles' reflection-induced force-derivative jumps do not invalidate this integral bound. For selected p take t=t_p and its incoming value; for nonselected p take t=t_g. Both satisfy `|s-t|<=Delta`, t<=b. Subtracting the actual kappa identity gives

```
B_p=|kappa_p|+Delta*|a0p+a1p|/2+b*D_a,p/2+J_p*b^2/4.
```

B_p bounds `|vminus_p-vphys,p(t_p^-)|` for selected p and `|vminus_p-vphys,p(t_g)|` for nonselected p. Reflecting the two incoming velocities preserves their difference. A selected particle accumulates at most A_p W more velocity change between its own hit and t_g. Therefore

```
xi_x,p=X_p
xi_v,p=B_p+A_p*W                       (selected)
xi_v,p=B_p                             (nonselected)
xi_a,p=B_R,p+L_p(xi_x)
field_error_g<=norm(rho,Inf)+M*sum_j |Qj|*xi_x,j.
```

These bound the original physical state at t_g against numerical x1/vplus/a1 at represented time s. There is no circularity: global constants, independent first-hit/prefix certificates and the AW/VW gates establish t_g; uniform X bounds force; B carries velocity. Remainder clearance has not been assumed.

## 5. Full original-trajectory remainder certificate

The true physical remaining time r=tau-t_g and numerical duration bn differ:

```
Rlo=tau-b,  Rhi=tau-a,   0<Rlo<=r,bn<=Rhi
|r-bn|=|s-t_g|<=Delta.
```

Compare physical flow from its actual state at t_g with the polynomial rooted at the **different** numerical state/time:

```
zR_p(t)=Q(x1p)+Q(vplus_p)*t+a1p*t^2/2
H_p(t)=xi_x,p+t*xi_v,p+t^2*xi_a,p/2+J_p*t^3/6.
```

Until p's next own hit, `|xphys,p(t_g+t)-zR_p(t)|<=H_p(t)` and its velocity difference from `vplus+a1*t` is <=`xi_v+t*xi_a+J*t^2/2`. These follow by integrating with nonzero initial errors. Let Hrem_p=H_p(Rhi).

For nonselected particles require both walls' quadratic inward-distance minima over [0,Rhi] to exceed Hrem_p. For selected particles require this for the opposite wall. Include endpoints and each interior stationary point in a quadratic minimum. A hypothetical first next hit within the actual r would contradict the strict margin and carried tube.

For a selected particle's same wall, let sigma=+1 left, -1 right. Its true distance at t_g is **nonnegative**, rather than reset to zero. Both terms below are valid lower bounds on its true inward speed, so their maximum is valid:

```
vL_p=max(nu_hit,p-A_p*W, sigma*vplus_p-xi_v,p)>0
accL_p=sigma*a1p-xi_a,p
F_p(t)=vL_p+accL_p*t/2-J_p*t^2/6.
```

Before any next own hit, `dphys(t)/t>=dphys(0)/t+F_p(t)>=F_p(t)`. Require F_p(0)>0 and F_p(Rhi)>0. Concavity makes F positive throughout and rules out same-wall return. Dropping the initial-distance term uses one-sided geometry; xi_x is still retained in the force/field uncertainty xi_a. A numerical inward velocity alone is insufficient.

Thus the initial certificates cover through the true cluster, and these carried tests exclude every additional/omitted physical contact through external time tau on the **same original external-start flow**. They supplement numerical-quadratic no-second-root checks, actual finite values and strict-interior final numerical candidates; physical clearance is not a substitute for numerical-stage validity. Failure or uncertainty rejects the whole candidate.

If A/V/J use horizon tau, apply the physical estimates only to actual t<=r, whose absolute time is <=tau. Requiring polynomial margins on larger [0,Rhi] is conservative and needs no physical continuation beyond tau. A separate claim about physical flow to t_g+Rhi would need enlarged-horizon constants, for example over tau+W.

## 6. Optional terminal bounds retaining time misalignment

Let actual stored final candidates be x2/v2/E2 after bn; a2 is its canonical exact gather. Measure

```
dx2_p=Q(x2p)-[Q(x1p)+bn*Q(vplus_p)+bn^2*a1p/2]
Z_p=max(|vplus_p+a1p*Rlo|, |vplus_p+a1p*Rhi|)
Xend_p=Hrem_p+Delta*Z_p+|dx2_p|.
```

Z bounds the affine polynomial velocity on the interval containing r and bn. The Delta*Z term is essential: it bounds zR(r)-zR(bn), not two states at an assumed common terminal time. Xend bounds the original physical x at external time tau versus x2.

Retain second canonical exact Istar2 and form the signed total defect

```
rhoTot=rho+[Q(E2)-Q(E1)+Istar2]
fieldEnd<=norm(rhoTot,Inf)+M*sum_j |Qj|*Xend_j
DaEnd_p=|alpha_p|*delta*norm(rhoTot,Inf)/h+L_p(Xend)
kappa2_p=Q(v2p)-Q(vplus_p)-bn*(a1p+a2p)/2
Vend_p=xi_v,p+|kappa2_p|+Delta*|a1p+a2p|/2
                    +Rhi*(xi_a,p+DaEnd_p)/2+J_p*Rhi^2/4.
```

Exact original-to-final endpoint kernel evaluation can equivalently provide rhoTot by telescoping. The invariant force decomposition and trapezoid proof establish these field/velocity bounds, with all stage/source and terminal-time errors retained. Norming separate signed rho defects first is merely looser; omitting one is invalid. Private admission bounds need not become new public/checkpoint arrays unless a later contract exposes them.

These are a posteriori external-step-start bounds, not a convergence proof. In particular J*t^2/4 is deliberately coarse; it does not itself establish a global second-order theorem. Isolated-impact accuracy still needs the smoothness/separation assumptions and independent remote refinement evidence above the precision floor.

## 7. Canonical Istar ownership and implementation limits

The [exact source kernel](../src/coupling/yee_charge.jl) remains the only
charge/path-geometry owner. Its private `_yee_charge_exact_sums` allocates the
same N+1 rational accumulator, invokes that kernel and returns the unrounded
values. `_yee_charge_store_exact_candidate!` performs the existing conversion,
finite check and candidate stores. The old exact-candidate wrapper still invokes
both in order and returns `nothing`; public validation, native selection and
publication stay in their existing owner. These helpers have source
implementation and essential test definitions; exact-revision runtime acceptance
is pending.

The reflection source candidate retains that same raw Istar in each local chunk
while using the canonical candidate conversion for its actual stored current.
Longitudinal impulse consumes the first N entries; the raw array has N+1
entries for the shared charge/current kernel. Float64 workspace sums and a rational lift of
rounded current cannot replace it. No A/g/deposition/path loop is copied.

All a0/a1/a2 come from the existing original-geometry gather/factor owners.
The [duration-free exact source row](../src/models/maxwell_vacuum.jl) and
[unrounded work value](../src/diagnostics/yee_diagnostics.jl) are implemented
private owners shared by the reflection source candidate. Its `_yi_carry!` and
`_yi_clear_remainder!` retain state, timing and signed field-invariant uncertainty;
they do not restart the physical proof from a zero-error numerical contact state.
Existing event-free methods retain their routing.

This contract presumes exact original operands, fixed valid sources, genuine first-hit certificates, physical nongrazing/ordering and AW/VW gates, s-in-hull, actual completed residuals, numerical remainder checks and every carried physical clearance test. No clipped interval, heuristic tie, exhausted refinement, unsupported second cluster or zero-error reset is success. It covers one external step from its stored original state, not a global entire-run initial-trajectory theorem or a point/sheath/contact/absorption model. Wall ledgers, transactions, restart, clocks and implementation resource limits
are described in the [model guide](YEE_IMPACT_LONGITUDINAL.md). Source enforcement
is implemented; exact-revision GitHub numerical, convergence, rollback/restart,
JET/compiler, formatting/CI and complete-step cost acceptance remain pending.
Absorption, multiple clusters per external step, transverse walls and backend
integration remain outside this candidate; CAP-002B remains open.
