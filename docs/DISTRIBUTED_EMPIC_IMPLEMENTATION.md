# Rank-owned periodic electromagnetic PIC

`distributed_empic(g,species; grid_dc=:project,...)` source-implements coupled
3D CPU Float32/64 EM PIC on actual PencilFFTs physical pencils. This first slice
is all-mobile, nonrelativistic CIC with one common particle/field step. It builds
on actual-pencil ownership/transport and the canonical vector Boris, midpoint,
current-projection and Maxwell expressions. It does not wrap a serial EMPIC or
replicate global fields/populations.

The implementation-first pass adds no tests or execution. Package load,
numerical/convergence, actual-rank rollback/liveness, compiler/JET, formatting/CI
and complete-operation memory/cost acceptance remain deferred. Matching-layout
complete checkpoint is a separately developed required slice B; until composed,
this core alone does not supply end-to-end restart or complete CAP-003.

## Construction, ownership and phase

Load PencilArrays/PencilFFTs with HybridPlasmaPIC and collectively construct a
`distributed_fourier_grid`. Every rank supplies an ordinary Vector with the same
ordered species count and q/m; local populations may be empty or have different
lengths. All particle leaves are independent ordinary Vectors with finite x/v,
nonnegative finite weights, finite signed q and positive finite m. Original
q/m factors feed the canonical mover. Nonrelativistic mode adds no |v|<c gate.
There is no mobile mask, shape/relativity selector, subcycling or species schedule.

The new extension-owned `DistributedEMPIC` carries E/B/J three-tuples and both
rho endpoints on its actual physical pencil. Seed owned E/B with
`PencilArrays.global_view`, then call `init_empic!(model,species)` with physical
full-time velocities. Init privately wraps/migrates, forms rho, corrects only
resolved longitudinal E, retains seed transverse/null E and seed B, zeros J,
copies the density endpoints and resets time/phase. B divergence is not projected.

Each particle belongs to the rank owning its wrapped CIC base node, using actual
Pencil ranges/topology and the existing binary-search partition lookup. Eight
CIC terms retain repeated nodes on tiny axes; only transport requests deduplicate
nodes. Remote stencil values/exact moment terms route directly to actual owners.
Migration can cross any number of domains and retains x/v/weight/id/tag/q/m.
Owner objects remain stable, but all component arrays can be replaced. References
to old leaves do not follow migration; IDs are not renumbered or certified unique.

A private per-species local count cache records the last published lengths.
Ready step, diagnostics and save reject changed local counts as well as changed
owner order/q/m. Cache values are local, not rank-equality configuration. Explicit
init may bind new counts and reset phase; migration updates the prepared cache.
Same-count replacement or source/velocity/weight edits cannot be detected and are
unsupported between ready operations. This is not a content/identity certificate.
There is no silent advanced-phase source-refresh API in this bounded EM slice.

Nearest rounding and gradual underflow are required. Converted dt is nonnegative
and either true stored zero or at least 2floatmin(T). A zero step validates and
returns without work/phase changes. Positive dt must advance the finite clock and
be strictly below the downward-rounded spectral `maxwell_cfl(model)`. The cached
bound is checked against its grid/speed derivation. It is a constant-step vacuum
spectral bound, not a plasma timestep prescription or variable-step stability law.

## Staggered common update

E and positions are full-time; after a positive step, velocity and B carry the
same negative preceding-half interval. First use primes velocity with the
existing physical acceleration backup and B by a backward magnetic half stage.
For preceding interval p and new h:

```text
Bc        = Bcarried - (p/2)*curl(Eold)
vnew      = canonical_Boris(vcarried, gather(Eold), gather(Bc), p/2+h/2)
xnew, mid = canonical_midpoint_drift(xold,vnew,h)
Jraw      = sum(q*w*vnew*S(mid))/cell_volume
J         = canonical longitudinal correction of Jraw and both rho endpoints
Bnew      = Bcarried - ((p+h)/2)*curl(Eold)
Enew      = Eold + h*(c²*curl(Bnew) - J)
```

These equations identify the existing native/range-completed local owners; they
are not a new synchronized vacuum KDK formula or exact-arithmetic trajectory.
The first prime and center remain distinct rounded stages. E gathers are reused
between first priming/kick when unchanged; B0 and Bc are separate samples. Each
species kicks, drifts and contributes current on every step. Positive phase
history is resolved and strictly sub-CFL; the unstepped boundary requires stored
zero J and identical rho endpoints. No time=step*h relation is imposed.

## Represented source map and limits

CIC charge/current use the existing exact moment callback: original q, weight,
velocity and represented stencil coefficients add on node owners before division
by the exact product of stored dx and one final T rounding. This is a deliberately
new distributed moment map, not serial normal-path density/current bit identity.
One scalar route owner handles charge and each current component; current uses
three routed passes, with the private midpoint borrow restored before communication.
No duplicate native deposition or full-population assembly is used.

Initial Gauss and current share a four-spectrum orchestration around the existing
modal projection. Current is corrected to the represented normalized endpoint
source; source preparation retains both original finite rho endpoints through
exact subtraction/division before its T rounding, rather than materializing an
overflowing T difference. All input spectra are consumed before any inverse/output
store. Common scaling has an upward-rounded finite conservative multiplier bound;
unsupported headroom, lost nonzero normalization, coefficients, modes or output
ranges reject. FFT/source rounding and general conditioning remain qualified.

The effective odd-Nyquist-zeroed wavevector is unchanged. Current keeps true DC
from raw deposition and zeros all three components at every other all-zero
DC/Nyquist combination. Mixed axis-Nyquist modes with another resolved axis remain
resolved/transverse under the canonical policy. Initial Gauss keeps seed null E;
no operator can enforce arbitrary charge outside its derivative range. rho is
retained rather than erased to make a full Gauss claim.

`neutrality=:require` uses original represented physical particle/background
charge with the declared tolerance. `:neutralize` declares its physical uniform
compensator. Required explicit `grid_dc=:project` separately chooses how represented
grid DC enters the constraint; standalone strict Poisson admission is unchanged.
Physical neutrality does not imply exact rounded grid neutrality. Diagnostics
separate physical residual/policy compensator, deposited-grid charge/DC compensator,
old/new grid-DC change, declared/stored volume and deposition discrepancies. These
are representations of one projection, not two added backgrounds. A grid-DC change
cannot be removed by a periodic divergence and remains visible.

Local magnetic/electric rows reuse canonical native order and exact range repair
on the finite represented curl. A nonfinite materialized curl rejects even if a
later timestep could rescue the completed field; the error identifies that boundary.
No original-field scaled-curl recovery, complete gather range closure, exact DFT,
cancellation-relative FFT, exact continuity/Gauss or discrete energy theorem is
claimed. Native gather/Boris/drift/local-field expressions retain their order for
equal inputs; the complete distributed map has different moment/FFT rounding.

## Collective publication, diagnostics and remaining work

Every rank calls the same collective APIs in the same order, including empty
populations. Grid metadata are read-only, grid/communicator remain caller-owned,
and concurrent use is unsupported. Storage must use exact CPU physical/spectral
pencils, expose overlap to Base.mightalias, and keep public/private fields,
particles, caches/midpoints, Fourier work and metadata independent.

Local validation, allocations, packet encoding/decoding, receive assembly,
source/field arithmetic and publication preparation agree before the next
communication. All physical fields, complete same-type replacement leaves,
count/scalar rows and phase/clock are ready before the first public store.
Commit performs fixed Array/Vector copies, typed leaf rebinding and Ref stores;
it does not resize, decode or revalidate. Recoverable earlier failure preserves
public state on every rank; scratch may change. MPI/transport loss, fatal resources
and storage that hides overlap retain the existing exclusions.

`electromagnetic_diagnostics(model,species)` reduces exact scalar charge/count,
kinetic, three momentum components and field energy, rounding totals once; field
energy is counted once with original c² and stored cell volume. It reports
nominal velocity/magnetic offsets -last_dt/2. Stored staggered energy is not a
synchronized or exactly conserved total. Absolute Gauss and continuity residuals
use represented filtered sources and materialized divergences, with source size
and unresolved density reported separately. Continuity retains DC and excludes
non-DC derivative-null combinations. The unresolved-density remainder also
contains filter roundoff. Unrepresentable scalar maxima/energies can return Inf;
no Inf/Inf relative residual is manufactured. At step0 continuity is `nothing`.
Diagnostics alter private work only and retain finite FFT/derivative boundaries.

The model owns 11 physical and 20 private real arrays plus one extra complex
spectrum, existing grid work, local staged particles/midpoints and temporary exact
moment packets. A normal positive step uses two curls and one four-spectrum
constraint: ten forward and nine inverse transforms, plus two charge and three
current routed moment passes, gathers and migration. First priming adds a B0
gather. Exact arithmetic, repeated routes, source normalization, derived-bound
checks and scalar reductions have unmeasured cost. No zero-allocation, speedup or
cluster scaling claim is made.

Complete restart must preserve particles, E/B/J/rho_n/rho_np1 and phase/clocks
verbatim through a distinct checksummed matching-layout schema. Count cache is
derived from admitted incoming leaves; stale destination values/counts may be
replaced without reading them as phase evidence. No field regeneration or migration
belongs in restore. Slice B is required before this family is end-to-end complete.
Relativity, fixed populations, other shapes, schedules/subcycling, collisions,
changed-layout restart, nonperiodic fields/particles, GPU and broad distributed
hybrid/completeness remain separate work.
