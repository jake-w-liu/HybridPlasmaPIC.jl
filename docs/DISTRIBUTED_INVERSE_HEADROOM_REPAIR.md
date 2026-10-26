# Distributed inverse-range regression

Pre-edit consumer brief (2026-10-27): `dealias_two_thirds!` on a finite 32³
sign-of-cutoff-kernel field returns Inf on the actual two-rank Pencil backend,
although an exact Rational Cauchy bound certifies every output representable.
The corresponding serial failure and both dependency inverse owners were
checked. The canonical full-operation headroom repair is in SpectralOperators;
HPIC adds an independent separable-convolution MPI regression and updates its
minimum dependency to the published repair. Existing global-volume/max scaling,
collective lost-nonzero checks, aliases, and staged publication remain the owners;
no transport or numerical multiplier is duplicated here.

Acceptance: Float32/64, shifted and uneven 32³/(32,24,16) data, min/latest Julia
and actual MPI1/2/4/8; independent complete-output representability and normwise
reference checks, existing operator/failure/resource suite and full compiler
comparison. Wrong local-volume, post-hoc Inf masking, double normalization and
serial-only fixes must fail. Source and existing first-filter repair remain
isolated in prior worktrees. Source evidence and final status are recorded under
`HybridPlasmaPIC-audit-evidence/2026-10-27/distributed-inverse-headroom`.

Local gate result: PASS on Julia1.12.7 and1.13.3. All MPI1/2/4/8 complete
operator suites plus the independent shifted/uneven convolution regression exit0.
SO passes7,738 assertions each. Ordinary144 profiles/432 samples per runtime,
including the transverse derivative caller, preserve exact0.1.4 full-state hashes
and0B warmed allocations. Full22 compiler blocks per runtime match0.1.4 exactly
after explicit source-root replacement only; Aqua HPIC9/SO10 pass. The final
local auditor verifies22 terminal0 commands, source pins and both actual two-rank
baseline failures. See `distributed-inverse-headroom/final-local-audit.json`.
The shared repair is registered as 0.1.5, tree
`12831f35e5026cded6879554f3001dc1b614d01b`, with a tagged GitHub release and
[all four dependency CI jobs green](https://github.com/jake-w-liu/SpectralOperators.jl/actions/runs/40374238781).
Root registry-only consumer integration also passes on both runtimes: only the
SpectralOperators Manifest entry changes from 0.1.4 to 0.1.5; actual MPI1/2
children each verify the active project, package root and registered tree and
pass all 28 new convolution assertions. Aqua and 24 concrete analyses pass per
runtime. Complete compiler aggregates and all 252/203 virtual-frame records
match the frozen parent after exact HPIC-root replacement only. Root evidence:
`takeover/inverse-consumer-integration/compiler-provenance-review-v1.json`.
Exact-SHA HPIC CI remains an acceptance gate.
No controlled timing, process-RSS or multi-node scaling result is claimed.

A separate invalid-alias robustness probe found that a same-pencil output can
be reshaped over a read-only real metadata vector on thin pencils and admitted
by the existing workspace validator. It has a separate planned validation slice;
this shared numerical headroom fix does not modify that admission behavior.

The full-operation bound applies when the common input scale is representable.
The existing `floatmax(T)` scale cap and rejection of nonzero inputs lost during
rescaling remain in force; a capped scale does not certify every representable
output when the required scale itself exceeds the format.
