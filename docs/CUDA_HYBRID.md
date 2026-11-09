# Resident CUDA isothermal hybrid

`cuda_hybrid_stepper` implements one periodic, single NVIDIA device, D1/3V
Float64 hybrid family. Ions retain all three velocity components; massless
isothermal electrons supply the quasineutral Ohm closure. Both transverse
magnetic polarizations and Hall evolution are included. This is a separate
concrete resident state under `CUDAHybridStepper`, selected explicitly after
importing CUDA. Existing CPU `HybridStepper` and Metal methods retain their
interfaces.

The source is implemented against the shared scalar interfaces and resident
Fourier primitives. Exact-revision GitHub and real NVIDIA acceptance must be
recorded separately. Source inspection does not establish successful kernel
lowering, device trajectories, memory use or speed. A CPU CI pass does not
certify the resident CUDA solver.

## Source and acceptance status

The complete family source has been reconciled against the current CPU scalar,
bit-admission, original-factor and checkpoint owners. SpectralOperators 0.1.9
contains the shared scalar Fourier products and D1 curl descriptor required by
this family. The CUDA executor deliberately retains its documented all-input
headroom policy; the serial operator's active-component selector is separate.

Actual NVIDIA execution remains required before accepting this solver. That
includes kernel lowering and separate RNE arithmetic, odd/even Fourier modes,
two transverse polarizations with Hall/pressure response, unequal-q/m populations
and an empty row, first/second steps with magnetic subcycling, late rejection
without public mutation, complete changed-count restart, and repeated-step device
storage/transfer measurements. No device trajectory, rollback, resource budget or
speed result is established by the source reconciliation.

## Verification entry points

`test/test_cuda_hybrid_host.jl` runs in a fresh hosted Julia process, separate
from the ordinary extension tests. The main CI matrix checks exact CUDA 6.4.1
on minimum Julia and supported CUDA 6.x on latest Julia. It imports the complete
extension and checks dispatch, dependency provenance and pure host admission;
it performs no resident evolution and establishes no NVIDIA numerical result.

`test/test_cuda_hybrid_device.jl` requires functional CUDA/cuFFT and has no
hardware skip. Its small odd/even cases compare complete resident evolution
with the CPU discrete method, exercise an empty species, complete changed-count
raw/archive continuation, and require a specific late private density rejection
with unchanged public state. These new test files are source-reviewed and remain
unexecuted at this integration checkpoint.

The manual `cuda-device-essentials.yml` workflow requires a registered Linux x64
NVIDIA runner carrying the `cuda` label. Run it on the exact candidate revision
and inspect all logs. The separate PTX/RNE, broader numerical/convergence and
steady storage/transfer/performance gates remain required; neither these small
fixtures nor a hosted import proves them.

## Load, initialize and advance

This API example is source-reviewed; no numerical result is asserted here.
CUDA6.4.1 or newer supported6.x, Julia1.12 or latest stable, functional CUDA and
cuFFT, and SpectralOperators0.1.9 or a supported later version are required.
The factory creates actual full-complex device plans and propagates failures.

```julia
using HybridPlasmaPIC
import CUDA

N, L = 64, 2π
np, n0 = 4N, 1.0
host = ParticleSet{1,Float64}(np; q=1.0, m=1.0)
load_lattice!(host, (0.0,), (L,), (np,))
fill!(host.weight, n0*L/np) # Number weight: canonical uniform-density normalization.

# Specify DeviceMemory explicitly; a CuArray default-memory preference can
# otherwise select unified/host memory, which this family does not admit.
dev(a::Vector{T}) where T = CUDA.CuArray{T,1,CUDA.DeviceMemory}(a)
ions = ParticleSet{1,Float64}(map(dev,host.x), map(dev,host.v),
    dev(host.weight), dev(host.id), dev(host.tag), host.q, host.m)
species = [ions]
B = (dev(ones(N)), dev(zeros(N)), dev(zeros(N)))
model = HybridModel(IsothermalElectrons(0.1); nfloor=1e-8)
st = cuda_hybrid_stepper((N,), (L,), model, species; dt=1e-4, NB=2, B=B)
init!(st)
step!(st)
snapshot = cuda_hybrid_host_snapshot(st)
budget = cuda_hybrid_host_diagnostics(st)
```

Select the fixed timestep/subcycle count to resolve the actual problem. The
example is not a stability or accuracy bound. The scalar particle convenience
`cuda_hybrid_stepper((N,), (L,), model, ions; ...)` delegates to the same
one-row collection convention. Positive steps must retain the saved dt;
`step!(st; dt=0.0)` validates the ready boundary and advances nothing.
Nonzero durations disappearing during Float64 conversion are errors.

## Model and numerical map

Units are Ωci-normalized, with normalized elementary charge and μ0 equal to
one. Each population retains its own original positive charge q, mass m and
number weights w. The grid density is the charge-equivalent electron density
`n=Σ q_s n_s`; ion bulk velocity is `Σ q_s n_s u_s/n`. Per-species CIC raw
moments are divided by stored dx, multiplied by q and combined in order before
one common normalization. Scatter and gather share the canonical CIC geometry.

The equations are `pe=Te*n`, `J=curl(B)`,
`E=-ui×B+(J×B)/n-grad(pe)/n`, and `Bdot=-curl(E)`.
Particles use the canonical nonrelativistic finite-step Boris/Euler rows.
The complete positive step primes once if needed, pushes/drifts all populations,
deposits frozen midpoint moments, performs all four RK4 stages in every magnetic
subcycle, deposits endpoint moments, computes a forward half-kick predictor,
and redeposits/recomputes the carried integer electric field. Predictor velocity
does not replace the carried half-time velocity.

The collocated Fourier operator uses the package's full mode metadata and
centered first derivative with even-grid Nyquist zero. Resident cuFFT transforms
are full-complex, normalized once on inverse, with canonical headroom scaling.
D1 curl includes all three input components in the common magnitude, including
the inactive longitudinal component. Constant Bx is retained without projection.
Atomic scatter order and cuFFT rounding need not match CPU bits or repeated runs.
No exact discrete energy theorem or general nonlinear stability is asserted.

## Admission, ownership and failure

This family requires CIC; finite Te≥0; positive normal nfloor; original
η=ηH=de2=0; `2≤N≤typemax(Int32)-2`; positive normal L/dx/dt and all required
half/subcycle/tableau durations. Original q,m and q/m must be positive normal
Float64. Physical Float64 inputs/results must be stored zero or normal;
unsupported native certificate, deposition replay, Ohm source-risk, RK4
completion or FFT rescaling-loss events reject. Actual density must remain
at least nfloor through every deposit; no floor activation or clipping occurs.

All particle leaves are ordinary, independent, contiguous DeviceMemory
CuVectors on one bound actual context/device, with Float64 x/v/weight, UInt64 id
and UInt32 tag. Empty populations are allowed within a nonempty ordered vector,
but initialization requires positive deposited density everywhere. Views,
wrappers, unified/host memory, cross-device leaves, aliasing and out-of-domain
positions reject. Kernel addressing follows explicit finite cell/index checks.

The state owns a copy of the ordered population-reference vector and execution
resources. Population objects/order and original q/m remain fixed. Counts and
loaded leaf bindings may change through explicit init/restore. Call `init!`
after loading new x/v/weights/B; it interprets velocities as physical v0 and
resets time/phase. Private configuration, metadata, work and plans must not be
modified externally. External writes during an operation are unsupported;
an atomic exclusive guard rejects overlapping public operations.

Successful mutations replace particle x/v and public field leaves. Refetch
`ions.x`, `ions.v` and `st.fields` afterward: old references are not snapshots
and may become private scratch. Ordinary steps retain weight/id/tag leaves;
restore replaces every particle leaf. Snapshots are detached ordinary host
arrays and remain usable. Do not externally resize/free/wrap shared empty raw
allocations while they belong to the state.

All expected validation/numerical failures occur on private candidates before
public particle/field/phase publication. Synchronization and caller context/
stream restoration precede the host pointer/phase commit. Private scratch may
change on failure. Fatal device/context loss, OOM and asynchronous driver errors
propagate; reuse of a poisoned context is not promised. The solver does not
destroy caller contexts or automatically reclaim the global memory pool.

## Phase, diagnostics and continuation

After init, x/B/E and loaded physical v are at time0, unprimed. After a positive
step, x/B/carried E are at integer time and live v is at nominal time-dt/2.
Advanced live/saved states require `dt≤time`, primed=true and last_dt=dt;
time is preserved as stored, without a time=step*dt requirement. Density,
velocity moments and pressure/gradient/reciprocal/current are the final carried
predictor fields and are part of the saved state.

`cuda_hybrid_host_snapshot` explicitly synchronizes and transfers every physical
particle/field/provenance leaf. `cuda_hybrid_host_diagnostics` takes one such
snapshot and invokes canonical host energy, momentum and Fourier divergence
owners with a separate host grid. Energy uses actual carried velocities and
integer fields/predictor density; velocity phase/offset accompany the result.
Neither operation silently reprimes, recenters or changes evolution state.

```julia
save_checkpoint("cuda-hybrid.chk", st)
archive_run("cuda-hybrid.run", st; rng_seed=0)
load_checkpoint!(st, "cuda-hybrid.chk")
load_archive!(st, "cuda-hybrid.run")
# Explicit species forms are also available; load_archive(path) only inspects.
```

Distinct schema1/policy1 `:cuda_hybrid_d1` persists all sixteen physical fields,
all ordered x/v/weight/id/tag/q/m records, full grid metadata, fixed model/dt/NB
and actual phase verbatim. The run archive uses the existing metadata/checksum
container; rng_seed is provenance, as this family owns no RNG.
Raw checkpoints use canonical atomic serialization. Only ordinary host payloads
are serialized; plans/streams/contexts/device pointers/private work are omitted.

Restore checks original links and strict policy/model/grid/phase/physical
admission, detaches/revalidates, prepares complete same-device replacements,
synchronizes, then publishes. Saved counts may differ; stale destination
physical values/phase are replaceable. No deposit, field solve, projection,
recenter or prime reconstructs saved values. Runtime/device versions are
provenance; compatible NVIDIA reconstruction does not imply bitwise continuation.
CPU/Metal/other-family payload migration is unsupported.

Hot init/step performs no full field/particle host transfer or host physics
solve. Bounded scalar status/max/min/certificate transfers are intentional and
synchronize the device. No zero-host-traffic, allocation or speed claim follows.
Real device/compiler/numerical/rollback/transfer and repeated-step memory gates
remain necessary before acceptance. D2/3, Float32, other shapes/closures,
variable stepping, full-range device cold repair, sources/collisions,
relativistic/ES/EM models, walls, Metal resident and multi-GPU remain separate
implementation tracks.
