# Case 27 — Hall-MHD oracle comparison with NHDS

This case compares the independent Hall-MHD eigenvalue oracle in
[`test/oracles/hybrid_dispersion_oracle.jl`](../../test/oracles/hybrid_dispersion_oracle.jl)
with output from the external NHDS kinetic dispersion solver. It does not measure
a wave evolved by the package's particle/field integrator. Production time-series
validation is tracked in [CAP-007](../../docs/PLASMA_COMPLETENESS_IMPLEMENTATION_PLAN.md#cap-007--independent-validation-and-comparative-evidence).

The reference uses NHDS's bundled `whistler.in` input. The local oracle selects
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

The builder requires Git, make and a Fortran compiler (`gfortran` by default;
set `FC` to override). The runner attempts a build if output is absent. Missing
output or too few parsed points produces a failed metric; malformed numeric
input can raise a parsing error. These failures are not skips. The case is
`default=false` and is outside the default package test suite.

Source, binaries and output are gitignored. The current builder clones upstream
HEAD and removes the clone's Git metadata; an existing executable is reused.
The runner reuses existing output without checking its source revision or input
checksum. Pinning, provenance verification and rebuild invalidation remain open
under CAP-007; current cached output does not establish a reproducible reference.
