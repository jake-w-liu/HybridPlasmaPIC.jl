# Case 27 — Hall-MHD oracle comparison with NHDS

This case compares the independent Hall-MHD eigenvalue oracle in
[`test/oracles/hybrid_dispersion_oracle.jl`](../../test/oracles/hybrid_dispersion_oracle.jl)
with output from the external NHDS kinetic dispersion solver. It does not measure
a wave evolved by the package's particle/field integrator. Production time-series
validation is tracked in [CAP-007](../../docs/PLASMA_COMPLETENESS_IMPLEMENTATION_PLAN.md#cap-007--independent-validation-and-comparative-evidence).

The reference uses NHDS's bundled `whistler.in` input, whose propagation angle
is 0.01 degrees. The local oracle selects
the highest parallel dispersion branch, with normalized background field and
density equal to one and sound speed `sqrt((5/3)*0.5)`.

## Acceptance and artifacts

[`case_27_nhds_dispersion_comparison`](run.jl) owns these gates:

- At least five points from output species/branch identifier `1`.
- Maximum relative frequency difference at `k*d_p >= 2` of 2%.
- Maximum relative difference over all parsed points of 8%.
- Maximum low-k difference at least twice the high-k maximum.

These are acceptance limits, not measurements from a fresh run. The last gate
checks an ordering of differences; it does not establish their physical cause.
The runner writes metric and frequency-overlay CSV files to its artifact directory.

## Running

```bash
bash validation/27_nhds_dispersion_comparison/build_nhds.sh
julia --project=validation validation/27_nhds_dispersion_comparison/run.jl
julia --project=validation validation/27_nhds_dispersion_comparison/plot.jl
```

The builder requires Julia, Git, make and a Fortran compiler (`gfortran` by
default; set `FC` to an executable path to override). `JULIA` selects the Julia
executable used by the shell wrapper. No tools are installed automatically.

Both commands use the shared [reference provenance owner](../reference_provenance.jl).
The [source lock](../reference_sources.toml) pins the full NHDS revision and
license checksum. Each attempt retains the detached source checkout, actual
input, build log, binary, raw output and a TOML checksum manifest. The runner
verifies this manifest before reading output and copies it into the validation
artifact directory as `27_nhds_dispersion_comparison_provenance.toml`.

New attempts live under `reference_cache/`; the older `nhds_build/` directory is
preserved and is not consumed. A changed input, builder or tool identity creates
a new attempt. Modified cached artifacts produce a visible failure. Preserve
such an attempt for inspection and move its cache directory aside to rebuild.
The shell builder also accepts an explicit cache-directory argument for isolated
reference generation. A writer lock prevents concurrent publication; a stale
lock must be inspected before removal. Failed attempts never replace the last
successful receipt.

Missing tools and failed build commands produce a failed availability metric
with the original diagnostic. Too few parsed points fail their existing metric;
malformed numeric input can raise a parsing error. These failures are not skips.
The case remains `default=false`, outside the default physics suite. Offline
provenance mutation/failure tests run in the package test suite.

The recorded identities support local cache verification; they do not constitute
a hermetic toolchain or guarantee bitwise agreement across machines. The physics
acceptance limits above are unchanged.
