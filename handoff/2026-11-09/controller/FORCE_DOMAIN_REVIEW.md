# Physical-v force and shape-domain review

Read-only derivation against immutable HybridPlasmaPIC123c4a23e59923f02ce96150c34ce32ccae72442 and inherited CONTRACT_V2. This is a bounded scalar-force/domain review, not production approval, numerical acceptance, or a full projected-pair proof. Only this evidence file is child-owned.

## Force derivation and canonical ownership

For physical p=m*gamma*v, dp/dt=q*(E+v×B), gamma=(1-|v|²/c²)^(-1/2):

    dp/dt = m*gamma*a + m*gamma³*v*(v·a)/c²,
    v·dp/dt = m*gamma³*(v·a) = q*(v·E),
    a = q/(m*gamma) * [E+v×B-v*(v·E)/c²].

These follow by differentiation and v·(v×B)=0. Q denotes exact represented-operand conversion, F=_boris_exact_factor(factor), and H=Q(h_actual), the ONE resolved finite positive T interval consumed by BOTH branches and their source/current projection. H is NOT the source-calendar gap delta unless equal; duration/public-clock defects remain independent receipts/tolerance gates.

For each stage j, define exact rational g_j=(Q(c)²-sum(Q(v_j)²))/Q(c)² and

    R_j = H*F*[Q(E_j)+Q(v_j)×Q(B_j)
                    -Q(v_j)*(Q(v_j)·Q(E_j))/Q(c)²],
    d_j = R_j*sqrt(g_j),
    v_low = v0+d0,     v_high = v0+(d0+d1)/2.

Stage1 uses the valid materialized Euler velocity, Euler full-time fields, its own gap and actual gathered samples. NR instead has R_j=H*F*(E_j+v_j×B_j), d_j=R_j. Independent symbolic checks: v=0 gives H*F*E; parallel E,v and B=0 gives H*F*E*g^(3/2); v·E=0 removes the radial subtraction. No numerical proxy was run.

VERIFIED source: _euler_kick_parts returns u +/- hq*(E+u×B), with parts4:6 the scaled kick (kinematics_range.jl:185–209). Calling it on all-Q operands reuses the canonical NR cross-product owner. It is not the relativistic physical-v RHS. _boris_rel_parts converts to momentum then applies the Boris map (full_pic.jl:699–734); _rel_prime_repair backs up momentum and divides by its new gamma (kinematics_range.jl:367–415). Substituting either finite map changes the requested embedded Euler/Heun method.

## Minimal proposed shared scalar API

These methods are proposed, not implemented:

    _em_lorentz_scaled_parts(v::NTuple{3,T},E::NTuple{3,T},
        B::NTuple{3,T},factor,H::Rational{BigInt},c::T,Val(rel))
        -> (coefficient::NTuple{3,Q},gap::Q)
    _em_lorentz_pair_output(v0::NTuple{3,T},p0,p1,c::T,Val(rel))
        -> NTuple{3,T}

p1=nothing requests Euler; second-stage parts requests Heun. A plain immutable parts packet suffices. The caller owns same-H coherence, gathered samples, coupled stages, transaction, clocks/calendars and budgets. Planar callers supply canonical absent E/B components as exact zeros and retain polarization admission (empic_collection.jl:851). Fixed rows do not force-step; zero-weight/zero-charge mobile markers still participate in domain/error observation.

- Scaled-parts checks finite original v/E/B/factor operands and positive finite c. Convert original operands to Q; obtain H*F*(E+v×B) via _euler_kick_parts parts4:6. NR returns that and gap1; rel subtracts exact H*F*v*(v·E)/c² and returns _exact_lorentz_gap, requiring gap>0 AND canonical supported input admission.
- Never form T(q/m), T(H*F), T(h/2), unscaled physical acceleration or separately materialized T increments. The ENTIRE NR Euler/Heun output is composed in Q, Heun divides by exact2, then _kinematics_round converts each completed component once. Specify preservation of original stored component bits when its exact net increment is zero, so a Q roundtrip does not erase structural signed zero.
- Rel composes v0+R0*sqrt(g0), or v0+(R0*sqrt(g0)+R1*sqrt(g1))/2, in ONE canonical _kinematics_bigfloat scope. Convert already-composed Q R/g through the exact conversion owner inside that scope; use _kinematics_round only on completed outputs. Do not round increments to T before combining. Exact-zero coefficient vectors can preserve v0 after input admission.
- Every materialized Euler predictor and high output must be finite and pass _supported_relativistic_output. A valid high endpoint cannot excuse an invalid force stage. Never clip/renormalize velocity or replace/widen C1 with exact gap positivity. A trial failure remains unpublished and retries/fails under bounded controller policy.

VERIFIED carrier: _BorisSpeciesFactors stores q,m,qm; _boris_exact_factor returns Q(q)/Q(m); _boris_factor_operands returns q,m (boris_accuracy.jl:4–16). _validated_boris_factor retains this IEEE CPU carrier even if derived qm is Inf or zero (boris.jl:439–445). Thus finite original factors must not be rejected merely because native qm loses range. An explicit scalar factor retains its represented-scalar meaning.

VERIFIED C1: rounded beta²<1 is required; only its admitted near-boundary side is resolved by exact stored squared norm (kinematics_range.jl:98–124; boris.jl:471–490). An exactly subluminal vector with rounded beta²>=1 remains excluded. The new helper must retain that exclusion at every stage/output.

Precision limitation: _kinematics_bigfloat is256-bit RoundNearest (kinematics_range.jl:29–39); _kinematics_round(BigFloat) rounds that already-approximated value via its exact represented rational (lines19–22). This is NOT a correct-rounding theorem for the exact radical expression, nor a general componentwise cancellation/midpoint, range-classification or true-error certificate. full_pic.jl:794 explicitly retains this fixed-replay caveat. Exact Q force/gap composition removes avoidable pre-rounding; independent high-precision/enclosure acceptance and any stronger certificate remain necessary. No accuracy/optimization measurement is claimed.

## Actual shape domains and strict-gate limit

VERIFIED coordinate: s_T=mod(x_T,L_T)/dx_T in caller T arithmetic (shapes.jl:83–101), not exact N*x/L. The latter bounded Yee CIC rule at line141 is a different owner. SpectralOperators-v0.1.9 stores dx=L/n (periodic_fourier.jl:83). Treat the physical periodic seam as a separate cut, including represented L/dx effects.

| Shape | Polynomial cuts | Canonical ownership |
| --- | --- | --- |
| NGP | k+1/2; value/source can jump | RoundNearestTiesUp: tie to k+1 (shapes.jl:116–120). |
| CIC | integer k; continuous value, derivative can jump | floor(Int,s), weights(1-f,f) (lines122–132). |
| TSC | k+1/2; value and first derivative continuous, second derivative can jump | round(Int,s), default ties-even; explicit quadratic weights (lines134–138). |

Regularity here describes real shape polynomials after periodic node identification; stored evaluation retains rounding. Exact-replay moments obtain the SAME stored canonical stencil, then lift its weights/factors to Q (_foreach_exact_particle_contribution, coupling/moments.jl:405–422; _exact_moment:432–441, read via git show at pinned commit). Complete gather likewise composes stored weights (gather.jl:63–90). A new gate must not assume a smoother exact x-to-stencil owner.

Concrete conservative event-free gate: retain actual represented lifted Euler/Heun positions and ideal/rounding receipts; bound travel/ambiguity and require one physical tile. For EVERY mobile marker/axis, require original/low/high physical positions divided by Q(dx) to lie in one open polynomial interval; additionally require actual canonical s_T/stencil selectors to match that piece after node mapping. Apply to all force/source stages and zero-weight markers. Each candidate straight segment lies in this convex interval, so it cannot exit/reenter unseen. Same wrapped endpoints or stencil bases alone do not prove this. Failed domain trials consume attempted budgets, restore the identical origin/calendar and halve only within explicit hmin/resolution bounds.

This is locally usable for strictly interior event-free trials, not a complete knot-crossing lifecycle or an exact-orbit certificate. Continuous CIC/TSC endpoint admission can later use a documented one-sided closed-piece/zero-extension identity, since footprints can change only by zero weights at the knot; mere selector equality is too restrictive there. It still needs explicit event/materialization/continuation acceptance for both branches and restart.

VERIFIED liveness counterexample by direct source/formula trace: NGP with exact dx=1, origin x=k+1/2 away from a seam, mobile neutral marker v_x<0, zero force. Ties-up selects k+1 at the origin. EVERY exact positive Euler displacement enters the left piece, so no positive shrink admits that physical trial under a strict open-piece rule. T rounding may stall at the knot; stalled represented motion is not successful passage. This proves only this admission limitation, not a general impossibility for NGP integration.

Full follow-through stays OPEN: CIC/TSC need approved event/right-edge/one-sided or separately justified crossing acceptance. NGP additionally needs a concrete jump-aware policy for actual changed endpoint charge, canonical continuity/DC/Nyquist projection, consumed field/current impulse, signed work, force side and complete event/restart state. A zero jump report or fictitious infinitesimal step cannot replace it. A deliberate new-policy shape restriction requires root approval and leaves full NGP error control open. A finite crossing can yield a measured low/high defect, but that alone supplies no smooth RHS order/effectivity theorem or halving-liveness guarantee.

## Acceptance and handoff

Force/API and bounded interior admission are ready for parent/root design review. General crossing, complete projected-pair consistency and scientific/effectivity/resource acceptance remain open. Essential independent REMOTE checks: NR exact-Q completed Euler/Heun with extreme original q/m and half-h underflow; relativistic enclosed/high-precision physical-v references including parallel/perpendicular fields and near-C1 inputs/outputs; fixed256 midpoint/cancellation witnesses; original-factor zero-field equilibria; exact/adjacent shape knots and physical seam, hmin/domain rollback, then approved crossing/continuation/restart. These checks are named, not run.

Only source reading, symbolic derivation, source-pin hashing and this evidence write were performed. No Julia/parser/compiler/formatter/numerical/proxy/test/device/profile workload, production edit, process control, CI, commit or push. Large caller coverage was bounded to cited owners; no whole-repository audit is implied.

## Checked source identities

Baseline: 123c4a23e59923f02ce96150c34ce32ccae72442. All supplied source manifest hashes and inherited contract hash match actual copies.

| Owner | SHA256 |
| --- | --- |
| src/particles/kinematics_range.jl | 5499b790c1691f1f5e48e9fdb9fb428fc7e6c3a1b32476c53a2513cc12df9df3 |
| src/particles/boris_accuracy.jl | 1432b8cab1c80b35bbe4d1b97d5110a9a6bee4c4853e564bfc685c451d6d03b6 |
| src/particles/boris.jl | 8403db9e8e8259bc353a783042a657053a64e51744e062195ea3a75d83d0a9e9 |
| src/models/full_pic.jl | 1b94500067c8baa28f3156dd9493aa49c8ab805c0e4d9e8f333f8327ac12fb60 |
| src/models/empic_collection.jl | 8e5abcdcffc116e66a7eabd0636972a94549502497b463760e906e74d4de9afd |
| src/coupling/shapes.jl | 39a65454ba01db1add461199fb4b7948c6878ebdda94a51384020298dfdede36 |
| src/coupling/gather.jl | f1b31e6dd7ad5ac538da27e4d99e2f82df2d5c8683e47e7d931d10b6156acd59 |
| src/utils/exact_arithmetic.jl | 9a49220b5d9286117114cd45ce896b70e0e1cb562f88c49127d5183531aa7b14 |
| git-show123c4a2:src/particles/boris_accuracy.jl | 1432b8cab1c80b35bbe4d1b97d5110a9a6bee4c4853e564bfc685c451d6d03b6 |
| git-show123c4a2:src/coupling/moments.jl | bbe261b3b5bafb961e94a53b9cba7babc5f469a5dbc8f17e38f5eadbcbcdf28a |
| SpectralOperators-v0.1.9/src/periodic_fourier.jl | 9d9efb2eb83574d28bdd1291013e6b1f135bc03b2667e6d524bb31778ffcbae4 |
