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

The builder requires Git, CMake and an MPI C/C++ toolchain. If `mpicxx` is absent
and Homebrew is available, the current script invokes `brew install open-mpi`.
It compiles the deck and launches one MPI rank. The runner attempts a build when
the profile is absent. Missing output or too few profile points produces a failed
metric; malformed numeric input can raise a parsing error. These failures are
not skips. The case is `default=false` and is outside the default package suite.

Source, binaries and output are gitignored. The builder clones the moving
`hybridVPIC` branch and removes its Git metadata. It reuses an existing compiler
binary and recreates its run directory. The runner reuses an existing profile
without verifying source, binary or input checksums. Pinning and provenance
verification remain open under CAP-007; current cached output does not establish
a reproducible reference.
