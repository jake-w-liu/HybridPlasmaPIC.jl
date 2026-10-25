# Validation frequency extraction repair

Status: local correctness, compiler, resource and case-20 checks pass on
Julia 1.12.7 and 1.13.2; push and exact-SHA CI remain pending.

## Pre-edit brief

- Mode: repair the canonical validation helper.
- Requirement: case 20's frequency magnitude comes from the dominant nonzero
  Fourier mode of its complex field history. A constant history has no measured
  oscillation and must fail validation rather than report an invented frequency.
- Reproduction: `first-bin-v1.jl` fails three of five assertions on Julia 1.12.7
  and 1.13.2. A 64-sample coherent first-bin tone at dt=0.25 returns
  9.8082943889 instead of 0.3926990817; zero and constant histories both return
  0.7853981634. Bins two and seven pass.
- Verified owner: `validation/common.jl:_peak_frequency` starts its search at
  bin two and never verifies nonzero spectral amplitude. Repository search
  finds one caller, `validation/20_empic_transverse_dispersion_2d/run.jl`.
- Compatibility: keep that caller and its 3% frequency / 1e-8 continuity gates.
  Preserve the existing short-history NaN result. Return frequency magnitude,
  including either rotation direction; record the exact estimator contract.
- Strategy: search the complete nonzero spectrum, map FFT indices to signed
  frequencies before taking magnitude, and reject constant/nonfinite signal
  evidence. Keep parabolic peak interpolation and explicitly remove DC from
  its adjacent-bin samples. No theory value enters the estimator.
- Interval contract: reject nonfinite/nonpositive `dt` with `ArgumentError`;
  the only caller supplies the finite positive constant 0.01. Unmeasurable
  histories return NaN, which the existing result owner records as failure.
  Center each component at its range midpoint and normalize the variation
  before the FFT; this keeps frequency independent of amplitude and DC offset
  without overflowing the sum of opposite finite endpoints. Test both ordinary
  and extreme amplitudes, including a tiny oscillation beside a large DC value.
- Verification: original witness, every coherent bin in odd/even lengths,
  both rotation signs, Nyquist, DC offsets, phase/amplitude variants, absent
  signal and off-bin accuracy against independently generated series; register
  tests in the package suite. Re-run actual case 20 on both supported runtimes,
  preserve raw histories/metrics, review formatting and exact-SHA CI.

Three inadequate repairs are rejected by these checks: changing only the first
search index leaves flat histories reporting frequencies; returning the expected
dispersion relation fails the synthetic-series matrix; searching only positive
bins cannot measure an oppositely rotating mode with the same physical frequency.
This peak estimator does not establish damping, uncertainty or multi-mode
separation. Those remain separate CAP-007B requirements.

## Progress

- [x] Read helper, caller, metric construction and documented validation scope.
- [x] Retain dual-runtime original failure and passing neighboring bins.
- [x] Verify both-sign/DC boundary witnesses before implementation.
- [x] Repair the owning helper and register independent regression tests.
- [x] Run actual production case, final source review and both runtime gates.
- [ ] Commit/push separately to main and review all six exact-SHA CI jobs.

Evidence: `../HybridPlasmaPIC-audit-evidence/2026-10-23/validation-frequency`,
including `first-bin-v1-commands.json` and complete logs.

Expanded baseline: `sign-dc-v2-commands.json` records six passes and 18 failures
on each runtime. The first driver had a wrong workspace path; its failed launch
and log remain retained and are not counted as a test result. The corrected
driver checks the source path before launch.

The signed-bin mapping follows the [FFTW transform definition](https://www.fftw.org/fftw3_doc/The-1d-Discrete-Fourier-Transform-_0028DFT_0029.html).
Julia's [official download page](https://julialang.org/downloads/) was checked
on 2026-10-25: latest stable is 1.13.2.

The amplitude-range baseline passes four and fails two checks on each runtime:
a separately represented tiny quadrature oscillation beside a maximum finite
DC component yields NaN (`baseline-range-v1-commands.json`). V1 centers the two
components before common scaling and passes all 3,153 registered assertions
per runtime (`functional-v1-commands.json`, source hashes in `source-v1.json`).
V2 reuses its private working array for the in-place FFT, removing the redundant
transform buffer while preserving the input history. `source-v2.json` pins the
exact final implementation, tests and caller; `final-v2-commands.json` retains
all six terminal commands and log hashes:

- 3,153 registered assertions pass per runtime; six formatter/AST/idempotence
  assertions pass.
- All four concrete compiler analyses have no findings per runtime. Forty
  numerical/input-preservation/allocation checks pass each. Every measured
  allocation count matches the parent for Float32/64 histories of length
  64, 256, 1,000 and 4,096; the first candidate used an extra array. These are
  allocation measurements, not latency or global optimization claims.
- Actual case 20 passes both unchanged gates on both runtimes: measured angular
  frequency 5.029486007282313 versus 5.0990195135927845 (relative error
  0.013636642520216225 < 0.03); charge residual 2.877175015333387e-12 < 1e-8.
  Raw histories and summary CSV files are retained in `case20-v2-{112,1132}/`.
- An independent input-type review initially hypothesized a Float16 regression.
  Its first probe failed in the unmodified parent: the installed FFTW does not
  support that input. The corrected comparison passes 12 checks per runtime:
  Float32/64, integer and rational inputs work in both owners; Float16 and
  BigFloat retain the parent's unsupported-type exception. No compatibility
  patch was made for that discarded hypothesis (`input-types-v2-commands.json`).

The canonical helper is the only behavior change; case 20 and its thresholds
are unchanged. Full package tests and complete compiler-tree review remain
required on the pushed exact revision. The broader CAP-007B campaign stays open.
