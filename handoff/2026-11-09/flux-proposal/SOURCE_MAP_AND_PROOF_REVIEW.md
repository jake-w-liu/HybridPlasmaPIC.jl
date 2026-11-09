# Source coverage and independent proof review

Read-only scope. Frozen 21 immutable base files plus four released CGL overlay files are recorded in SOURCE_PINS; a pin is not a claim of full-file audit. Actual relevant reading in this task:

| Source | Read coverage / purpose |
|---|---|
| electrons/closures.jl | 1–340, relevant scalar types/traits/grid conversions; current CGL overlay remains separately frozen, root-owned |
| electrons/energy_equation.jl | Entire 596-line file: pressure/heating, complete sum, quotient, dyadic preparation, full-k source and storage boundary |
| electrons/field_aligned_conduction.jl | Entire 162-line operator: workspace, original projector, channel preparation/mean projection, both outgoing channels and adapter |
| integrators/extrapolated_leapfrog.jl | 1–455,540–705,885–1241,1390–1605,1660–1745: allocation/alias, transaction, RHS/RK, current/force/phase, publication; broader energy routines were not freshly audited in full |
| models/hybrid_pic.jl | 1–110, transport/model/density/field boundary; later deposition was not freshly re-audited in full |
| diagnostics/conservation.jl | 365–440, fresh charge-density and instantaneous conduction contract; prior thermal diagnostics are preserved |
| io/checkpoint.jl | 230–335,375–675,775–855,1000–1157,1465–1498: field/model policy, complete graph, capture, actual restore and archive delegate |
| io/collisional_checkpoint.jl | Entire 339-line base; released CGL delta separately pinned/read previously, no new source change |
| integrators/collisional_hybrid.jl | 65–160,320–585 plus referenced owner inventory; released CGL delta separately pinned/read previously |
| utils/exact_arithmetic.jl | 1–80,309–399, exact extraction and directed root; no runtime/parser validation |
| electrons/thermal_exchange.jl | 30–140,165–230, root/old config/mass binding; root already owns integrated A/B source |
| verification/normalization.jl | 1–75,110–205,348–393, actual pressure/velocity/mass scales and reference transverse dimensions |
| HybridPlasmaPIC.jl | 70–175 include order; source-wide actual trait/caller grep for supported/rejected families |
| Documentation | Full current field-aligned operator policy considered; completeness plan 991–1041 and1043–1174 checked; existing CGL guide95–143 checked |

Other pinned force/state/integrator files provide immutable caller context, not fresh complete audit coverage. Source-wide searches located no current flux-limited closure/owner at this immutable base. No claim that arbitrary code outside this read scope is bug-free.

Settled algebra checked against these actual owners:

1. Normalization: actual `vA=B0/sqrt(mu0*n0*mi)`, pressure scale `B0^2/mu0`, and mass scale mi imply `p0=n0*mi*vA^2`. CM77's electron flux divided by p0*vA gives the contract's alpha*pe*sqrt(Theta/mu_e). Omitted coordinates keep their reference transverse measure; no new collision macro-count factor enters.
2. Channel adjoints: ideal real odd D is skew-adjoint and its stored Nyquist row is zero; N is self-adjoint nonnegative. The code independently applies both. For q=-aG, H=-Dq_D+Nq_N, pairing with Theta gives sum G*q=-sum a*G^2. Constant kappa and a->kappa gives D^2-N^2=-kfull^2 modewise. These are ideal metadata/operator identities, not a measured FFT theorem.
3. Actual line-mean projection: each computed outgoing channel has exact Q line sum zero by its literal total/length subtraction. Hence the sum of computed H is exactly zero before complete native publication, regardless of finite FFT roundoff. State/rate publication and RK addition are separate rounding operations.
4. Coefficient enclosure: kL<=sqrt(K)<=kH and rL<=sqrt(K*g2/C)<=rH imply aL<=a<=aH. Then aL^2*g2<=K*g2/(1+sqrt(K*g2/C))^2<=C. All zero branches precede division. This certifies constitutive rational q, not its later dyadic preparation.
5. Required preparation distinction: channel input is dyadic2(value-anchor), with unchanged `_conduction_dyadic` residual bound u^2*|value-anchor|. Restoring each axis's discarded constant anchor gives actual represented flux q_tilde. A local cap may be arbitrarily smaller than another point's anchor; the constitutive cap is therefore insufficient to prove q_tilde's cap from this error bound. The proposed strict exact prepared-input cap check is an unapproved admission choice, not a present-source defect or measured failure.
6. Entropy: for fixed actual n, differentiating n*log(pe/n^gamma)/g with respect to pe gives 1/(g*Theta). The actual source gate is V*s/g sum(H_computed/Theta), not variance. Pi + chain defect + measured adjoint defect cancels exactly to that sum, including projection/preparation. The rowwise reciprocal chain expressions have different dependence on an independently varying Theta_j, so no generic continuum chain identity has been established for this collocation map.
7. Publication: Q(out)=R_other+sH+delta by definition at the canonical complete-RHS owner. Keeping all three terms makes the entropy/energy rounding contribution explicit; a source-only gate proves no full RK, changing-density or plasma entropy/conservation theorem. No negative-entropy runtime fixture, error frequency or performance result was measured.
8. Frozen-mobility timestep estimate: with fixed positive n and a bounded by kappa_max, the ideal operator on p is s*L_a*diag(1/n). It is similar to a symmetric negative-semidefinite matrix through diag(sqrt(n)); its spectral radius is at most s*kappa_max*sum(kfull_max^2)/nmin. This justifies only a frozen-coefficient RK4 estimate. Derivatives of state-dependent a and the combined coupled map are outside that proof.

Unexecuted gates and decisions are listed explicitly in PRE_EDIT_CONTRACT. Literature coverage limitations are in PRIMARY_SOURCES. Root's stop request is honored: these are evidence files only; no production/test/shared docs, local numerical workloads, CI, commits/pushes or process controls were performed.
