# Distributed operator takeover audit (2026-10-27)

Baseline: HPIC `cd558315956df0debb14b4a2fb7790eaef2821e6`; registry-resolved
SpectralOperators 0.1.3, tree `036881a8627ae86e285825891bb427dd66508c08`.
Latest stable Julia checked at the official downloads page: 1.13.3.

## Pre-edit repair brief

The accepted infinite-strength filter produces NaN for positive mode powers
that underflow (reproduced in serial and distributed Float32/64 at p=10000).
The canonical repair belongs to shared SpectralOperators `_apply_exp_filter_multiplier!`;
HPIC adds only a consumer regression and a dependency lower bound after publication.
The shared-kernel repair brief is `SpectralOperators/docs/EXPONENTIAL_FILTER_RANGE_REPAIR.md`
in its isolated repair worktree. The intended result is the preserved DC mean,
with all other modes removed, even at the maximum accepted power. Existing
finite-strength behavior, rank agreement and rollback requirements stay intact.

Wrong repairs rejected by the checks: dropping the mean, restricting the already
accepted infinite strength, fixing only one large power, adding a consumer-only
formula, or accepting only the one-rank result. No production orchestration has
been changed in this consumer patch.

## Audit matrix

| Surface | Checked owner/invariant | Status |
| --- | --- | --- |
| Dependency/runtime | Registry source tree, minimum/latest Julia | 0.1.4 registered tree `1286d89fc8b05aa95639d6be5c432208f7296183` verified; registry-only 1.12.7/1.13.3 resolution passed |
| Logical indices/layout | Physical copies share storage order; mode kernels consume global logical views | Source traced and MPI1/2/4/8 analytic/serial/uneven-grid suites passed |
| Scaling | Global maximum over ranks/components, global transform volume, canonical gain | Boundary/global-scale suites passed; independent inverse-overflow defect repaired and registered in 0.1.5; its consumer regression is a separate integration |
| Aliases/publication | Exact physical pencil, CPU backing, private workspace disjointness; all components staged | Existing MPI1/2/4/8 vector/publication failures and later liveness passed; geometry-metadata alias admission is a separate confirmed repair in progress |
| Numerical operators | Shared derivative/curl/Laplacian/projection/filter kernels | Infinite-strength defect repaired in registered 0.1.4; inverse-range consumer integration and metadata alias admission remain separate |
| Resources/compiler | Reused buffers, direct-transform allocation control; full compiler trees | MPI allocation controls passed;128 shared-kernel profiles/runtime have baseline-identical hashes/0B;22 full compiler blocks/runtime match baseline |
| Bounded Maxwell integration | Vacuum and prescribed 1D current remain separate CPU Yee storage | Docs and call chain read; particle coupling, MPI/GPU and restart explicitly open |

Evidence and all failed harness attempts are retained under
`HybridPlasmaPIC-audit-evidence/2026-10-27/distributed-takeover`.

## Checked normal-path communication and storage inventory

`_mpi_collective_stage` owns exactly one success-path Allreduce. Tracing the
actual orchestration gives the following explicit scalar reductions (excluding
the transpose traffic inside PencilFFTs); the counts are source-derived, not a
multi-node performance measurement.

| Operator | Forward / inverse transforms | Allreduces | Broadcasts |
| --- | --- | --- | --- |
| derivative | 1 / 1 | 6 | 1 |
| gradient | 1 / 3 | 9 | 0 |
| divergence | 3 / 1 | 9 | 0 |
| curl | 3 / 3 | 11 | 0 |
| Laplacian | 1 / 1 | 5 | 0 |
| projection | 3 / 3 | 9 | 0 |
| exponential filter | 1 / 1 | 6 | 2 |
| two-thirds filter | 1 / 1 | 5 | 0 |

The grid explicitly owns one physical complex array, three spectral complex
arrays, and three physical real publication arrays. Their combined payload is
`sizeof(T)*(5*local_physical_cells + 6*local_spectral_cells)` per rank; summed
over ranks it is `11*sizeof(T)*prod(n)`. Physical and spectral local volumes can
differ under uneven partitions. Replicated axis vectors occupy
`(sizeof(Int)+4*sizeof(T))*sum(n)` payload bytes per rank. These formulas exclude
the dependency-owned FFT/transpose buffers and object headers; they are not a
claim about complete process RSS. The plan owns its temporary input/output
buffers and references its pencils/topology for the same lifetime. A caller's
communicator is retained and no HPIC operator frees it.

The checked `PencilArrays.global_view` wraps the original PencilArray in an
OffsetArray, uses global logical ranges and corrects permutation. It does not
copy the field. Its `Base.dataids` delegates to actual backing storage, which
makes the private-workspace and component alias checks meaningful. The consumer
extension contains no field-gather call in its operator implementation.

## Bounded first-repair verification

The canonical infinite-strength repair and consumer regression pass18 terminal
commands across Julia1.12.7 and1.13.3: full SO suite7,674 assertions each;
all8 MPI1/2/4/8 operator launchers;128 finite-strength profiles and384 samples
per runtime with exact baseline hashes and0B; full22 compiler blocks per runtime
exactly match baseline after explicit source-root normalization only, including
frames/signatures/messages. Aqua HPIC9/SO10 and formatter AST/idempotence checks
pass. Evidence: `distributed-takeover/final-local-audit.json` and `HANDOFF.md`.
All four dependency CI jobs pass at `0e17596`; 0.1.4 is registered and has a tagged GitHub release. Root registry-only resolution selects exactly the tested tree on Julia 1.12.7 and 1.13.3; all other dependency entries match the pre-integration environment. Root integration passes all four registry-only MPI 1/2 launchers on Julia 1.12.7/1.13.3. Both child ranks independently verify the loaded HPIC root, project, dependency version/tree and registry tracking. Aqua and all nine concrete analyses pass per runtime; the complete 29/27-report aggregate compiler forests, including 197/152 rendered frames, are byte-identical to the accepted parent without normalization. Exact-SHA HPIC CI remains the publication gate.

Independent audit then reproduced another defect on a finite32³ cutoff-kernel field: the
unnormalized inverse overflows before its final1/N normalization, although exact Cauchy bounds
certify every output representable. Both serial and actual two-rank backends reproduce it. This has
a separate shared-headroom repair and version 0.1.5 campaign. That dependency repair has all four CI
jobs green, is registered with tree `12831f35e5026cded6879554f3001dc1b614d01b`, and has a tagged
release. Its HPIC consumer regression will be integrated separately. CAP003 is not complete, and no
controlled latency or multi-node scaling result is claimed.