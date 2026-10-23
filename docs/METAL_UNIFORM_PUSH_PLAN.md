# Metal uniform-field particle push

Status: 2026-10-23; parent stack overflow reproduced on an actual Metal device.

## Pre-edit brief

- Mode: repair the Metal uniform-push dispatch path.
- Requirement: advance every device particle using the canonical Boris velocity
  update and then drift its D spatial coordinates; return the original particle
  set. Validate host parameters before any device writes, including empty sets.
- Baseline: `metal-uniform-push/baseline-v2-1132.log` records a valid initialized
  particle, a passing CPU free-drift control, and `StackOverflowError` on actual
  Metal 1.11.1 / Julia 1.13.2. Relevant parent source blobs are checked against
  `7b50b51`. The initialized minimum-Julia witness reproduces the same failure.
- **VERIFIED cause:** the extension's `push_uniform!` at line570 calls
  `_push_uniform_broadcast!`; `src/particles/boris.jl:210` calls `push_uniform!`
  again and selects the same extension method. The runtime trace repeats this
  cycle 39,991 times. This is independent of the separate NGP shape repair.
- Contract: preserve the public signature, Float32 backend restriction, particle
  identity, weights/tags/IDs, kick-before-drift ordering, finite host parameter
  validation and asynchronous Metal launch behavior. No CPU arithmetic changes.
- Canonical owners: reuse `boris_kick` and the existing CPU charge/mass, timestep
  and uniform-field validators. The Metal extension owns thread dispatch and
  device indexing. Remove the private recursive helper after its sole caller
  migrates; repository-wide search found no other caller.
- Compatibility: the broken device path gains real behavior. The private helper
  is not exported. CPU calls keep their existing scalar owner and semantics.
  Retain the extension's Int32 particle-index limit with explicit host rejection.
- Strategy: one in-place device kernel, one thread per particle, with no temporary
  device arrays. Check the unsigned global thread index before converting to
  Int32, including the padded last thread group. Reuse `_metal_launch_config`.
- Verification: actual Metal on both supported Julia runtimes; initialized
  free drift and electric-only oracles, gyro invariants, full CPU trajectory
  comparisons, 0/1/255/256/257/4096 particle counts and D1–3; invalid-input state
  preservation; full GPU/particle neighbors, compiler and formatting checks.
  Observe warmed host allocation and synchronized elapsed cost without claiming
  a speedup over a parent path that does not terminate normally.

## Progress

- [x] Reproduce the valid latest-Julia parent failure and trace the dispatch cycle.
- [x] Complete minimum-Julia reproduction (one CPU control pass, one device error).
- [x] Implement the canonical device path and regression tests.
- [x] Pass actual device, neighboring, validation, compiler and format gates.
- [x] Review measured resources and the final source/held-file hashes.
- [x] Commit and push this repair separately as `bf2a8bc`.
- [ ] Require all six exact-SHA CI jobs.

Wrong patches to reject: invoking the CPU particle loop would require prohibited
scalar device indexing; updating velocity components in separate dependent
kernels can read already-modified components; omitting the padded-thread guard
accesses beyond the final particle. Scalar-index-disabled execution, independent
gyro/electric oracles, and thread-group boundary sizes distinguish these cases.

The first ad hoc witness left its initial position uninitialized and continued
checking after the overflow; its failed logs remain retained but are not the
valid acceptance baseline. V2 initializes position and aborts the testset on an
exception before any post-overflow device access.

Evidence: `../HybridPlasmaPIC-audit-evidence/2026-10-23/metal-uniform-push`.
The NGP integration failure remains retained in `ngp-nearest-node-combined`.

## Local acceptance

Both Julia 1.12.7 and 1.13.2 use actual Metal 1.11.1 with a required functional
device and scalar indexing disabled. Each passes 778 focused assertions and
894 complete GPU-backend assertions in fresh processes. Coverage includes
all three spatial dimensions, empty/partial/full thread groups, independent
free/electric and gyro oracles, ordinary CPU trajectories and rejection before
mutation. The valid minimum/latest parent reproductions remain failed as
expected; no checks access device state after their stack overflow.

Each CPU lane passes 57,583 neighboring particle/validation assertions and
14 concrete compiler checks. Complete package JET trees match the accepted
parent after verified source-line mappings (30 / 27 existing reports); no
report was suppressed or classified as harmless. Formatting passes eight
syntax-preservation/idempotence assertions. The final verification checks
274 frozen Julia files, 355 unchanged committed Julia files, all 14 command
logs/driver hashes and the held radiation hashes.

Each runtime also passes 60 state/storage checks over 15 resource profiles:
D1–3 and 0/1/257/4096/262144 particles, nine warmed samples each. Empty pushes
allocate 0 host bytes; nonempty samples allocate 1,072–1,488 bytes on minimum
Julia and 1,056–1,472 bytes on latest. Every particle array retains identity.
Per-profile median synchronized times span 0.399–0.882 ms and 0.505–1.099 ms,
respectively. These observations were collected alongside untimed CPU work;
they are neither isolated latency acceptance nor a speedup comparison against
the failing parent. Device memory-pool counters are unavailable through the
current backend telemetry; no measured device-allocation bound is claimed.
The kernel source has no temporary array construction and launches once for
each nonempty call.

Evidence files: `final-verification-v1.json`, `device-commands-v1.json`,
`commands-112.json`, `commands-1132.json`, `resource-commands-v1.json`, raw
resource TOML, and `compiler-comparison.json` under the evidence path above.
GitHub CI remains required after the separate push; its Ubuntu jobs do not
replace the local actual-Metal execution.

Pushed as `bf2a8bc`; all 359 committed Julia files match the validated candidate.
[Exact-SHA CI](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/39909168466) is running.
