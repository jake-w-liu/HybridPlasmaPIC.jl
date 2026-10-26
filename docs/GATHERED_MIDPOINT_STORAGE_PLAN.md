# Gathered-push midpoint storage

## Repair contract

`push_gathered!(ps, E, B, dt; xmid)` publishes final particle coordinates and
optional midpoint coordinates. Midpoint outputs must be writable vectors,
disjoint from particle storage, gathered input fields and one another under
the canonical particle-storage alias predicate. Both CPU and Metal must reject
a conflict before changing particles or launching a device update. Sharing
read-only E/B components remains allowed. `xmid = nothing` retains its existing
path. This repair adds no buffers, particle fields, archive state or dependency.

**Mode:** repair. The retained one-particle witness starts with x=1, v=2,
zero fields and dt=0.5, passing `xmid=ps.x`. CPU and actual Metal both produce
final x=2.5 instead of 2, because the midpoint write changes the coordinate
used by the final drift. This violates the two-output contract even though
both outputs cannot coexist in the same cells. Rejection is the compatibility
decision for such overlapping buffers; no copied-field fallback is introduced.

**Owner and affected surfaces:** `_check_gathered_field_storage` in
`src/particles/boris.jl` already precedes mutation in both the public CPU
method and `ext/HybridPlasmaPICMetalExt.jl`. The storage preflight
check uses `_mightalias_particle_storage` and `_particle_arrays_mightalias`
from `src/particles/particle_set.jl`. The CPU native/guarded arithmetic and
Metal broadcast kernels remain unchanged. Current repository callers are the
public methods and their tests; coupled solvers use their own canonical
push/gather paths.

`Base.mightalias` detects overlapping CPU views, including distinct view
objects. The particle predicate follows parent-linked CPU wrappers, including
empty ones. Actual device behavior is a separate acceptance gate: a measured
nonempty contiguous Metal wrapper is detected, while a distinct empty native
Metal wrapper sharing the same device allocation is not. An empty wrapper has
no particle elements to overwrite. Detection does not cover every custom array
or unsafe independently wrapped pointer. Read-only midpoint storage remains
outside the writable-output contract.

## Acceptance tracker

- [x] Preserve the original CPU and actual-Metal witness and pre-edit sources.
- [x] Canonical preflight rejects conflicts with every particle storage family,
  each E/B component, other midpoint components and overlapping views.
- [x] Rejections preserve particle state, gathered fields and midpoint buffers;
  include empty and zero-step cases and validate the reported argument.
- [x] Valid disjoint buffers, disjoint CPU views and shared read-only E/B inputs
  retain their output and allocation behavior against frozen parent controls.
- [x] Formatter AST/idempotence, concrete analysis, Aqua and full compiler trees
  on Julia 1.12.7 and 1.13.3; actual Metal parity and preflight rejection.
- [x] Root integration repeats CPU, actual-Metal, resource and compiler gates.
- [x] Separate main push as `a79a07b`.
- [ ] [All six exact-SHA CI jobs](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/40400955427).

Each runtime passes 5,230 CPU and 416 actual-Metal assertions. All 54 matched
public state/allocation profiles retain their parent results; measured byte
counts are 0, 16 or 32 depending on the profile. Eighteen concrete calls and
Aqua pass. Full package analysis retains the parent's 29/27 reports and
252/203 virtual frames. The isolated comparison accounts for four exact
generated-function label occurrences in byte-identical owning files; root
integration matches that isolated candidate with only the checked source-root
mapping. Numerical arithmetic outside the admission function is AST-equivalent
to the parent after the helper rename.

The first root Metal setup changed an unrelated Krylov version and was held
before numerical tests. The corrected setup preserves all dependencies while
selecting registered SpectralOperators 0.1.5. Both failed setup and successful
ten-command integration receipts remain in the external evidence archive.

Checks must reject identity-only admission, a positions-only alias check,
validation after velocity publication, and a CPU-only fix. This is an alias
admission repair, not a universal transaction guarantee for arbitrary custom
arrays or an allocation/latency optimization claim.

Evidence is retained under the dated external audit archive at
`takeover/gathered-midpoint-alias`; historical primary work remains untouched.
