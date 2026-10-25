# Case 28 — Perpendicular-shock comparison with Hybrid-VPIC

This case runs the committed [`shock_perp.cxx`](shock_perp.cxx) deck with external
Hybrid-VPIC, extracts scalar shock metrics from its magnetic-field profile, and
compares them with `run_perp_shock_rh` at the inferred shock-frame Mach number.
The local run uses an RH-initialized shock. This comparison covers compression
and overshoot; a shock-aligned profile comparison and broader model validation
remain in [CAP-007](../../docs/PLASMA_COMPLETENESS_IMPLEMENTATION_PLAN.md#cap-007--independent-validation-and-comparative-evidence).

## Acceptance and artifacts

[`case_28_hybridvpic_perp_shock`](run.jl) owns these gates:

- At least 20 external profile points.
- External overshoot within `[1.1, 1.5]`.
- External compression no greater than the computed fluid RH compression.
- Relative overshoot difference between the two runs no greater than 25%.
- Relative compression difference no greater than 15%.
- Inferred external shock-frame Mach number within 5% of 6.

These are acceptance limits, not measurements from a fresh run. Passing these
scalar checks alone does not validate the complete plasma model or explain a
difference between the simulations. The runner writes metric and external-profile
CSV files; the plot displays the external profile.

## Running

```bash
bash validation/28_hybridvpic_perp_shock/build_vpic.sh
julia --project=validation validation/28_hybridvpic_perp_shock/run.jl
julia --project=validation validation/28_hybridvpic_perp_shock/plot.jl
```

The builder requires Julia, Git, CMake and a working MPI C/C++ toolchain.
`MPICC`, `MPICXX` and `MPIEXEC` select executable paths; `JULIA` selects the Julia
executable used by the shell wrapper. No tools are installed automatically.
The committed deck is built through upstream CMake's `USER_DECKS` hook and run
on one MPI rank. Release flags include `-fno-strict-aliasing`, as advised in the
pinned upstream README. The source tree remains unmodified.

Both commands use the shared [reference provenance owner](../reference_provenance.jl).
The [source lock](../reference_sources.toml) pins the full Hybrid-VPIC revision
and the checksum of its README, which contains the BSD-3-Clause license.
Each attempt retains the source checkout, deck and `injection.cxx`, actual build
commands/log, compiled program, raw outputs and a TOML checksum manifest. The
runner verifies the manifest before reading the profile and copies it into the
artifact directory as `28_hybridvpic_perp_shock_provenance.toml`.

New attempts live under `reference_cache/`; the older `vpic_build/` directory is
preserved and is not consumed. Changed input, builder or tool identities create
a new attempt. Modified cached artifacts fail visibly. Preserve such an attempt
for inspection and move its cache directory aside to rebuild. The shell builder
also accepts an explicit cache-directory argument. A writer lock prevents
concurrent publication; inspect a stale lock before removal. Failed attempts
never replace the last successful receipt.

Missing tools and failed build commands produce a failed availability metric
with the original diagnostic. Too few profile points fail their existing metric;
malformed numeric input can raise a parsing error. These failures are not skips.
The case remains `default=false`, outside the default physics suite. Offline
provenance mutation/failure tests run in the package test suite.

The recorded identities support local cache verification; they do not constitute
a hermetic toolchain or guarantee bitwise agreement across machines. The physics
acceptance limits above are unchanged.
