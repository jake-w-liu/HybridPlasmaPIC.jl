# Parallel implementation workflow

Updated 2026-11-09. **Codex stops at the user's requested OpenCode handoff. All workers have finished.** The full simulator goal remains incomplete; OpenCode has not started. No exclusive writer assignment carries over implicitly.

The [OpenCode plan](OPENCODE_IMPLEMENTATION_HANDOFF.md) owns P00–P43, dependency ordering, parallel lanes and G1–G10 acceptance. The [79-row tracker](OPENCODE_ACCEPTANCE_TRACKER.md) maps every original requirement; the [original completeness plan](PLASMA_COMPLETENESS_IMPLEMENTATION_PLAN.md) retains milestone history. The [audit tracker](DEEP_DEBUG_TRACKER.md#current-progress) retains checked revisions and historical failures. The [packet index](../handoff/2026-11-09/README.md) identifies all frozen work and hashes.

## Starting revision and blockers

Source baseline is `5aa628513fb854c96b2b387780edef1ea1faa66b` on `implement/cgl-integration-20261109`; the final handoff adds docs and evidence only. Accepted-main snapshot is `9b4d363806722319c2241e4af04ed1ed1da00119`, whose [26-job run](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/41813949753) is green. The unaccepted source stack has not been promoted to main.

P00 verifies frozen artifacts, primary dirty-work hashes and the compiler ledger in an isolated checkout. P01 repairs the verified parser failures in `src/models/empic_collection.jl:1430/1854` first, then classifies all complete added JET trees and fixes actual format/focused/MPI failures. The [priority register](../handoff/2026-11-09/ISSUE_REGISTER.json) names the checked evidence. No actual billing exception is established.

CGL source and its 238-line essential are integrated, with new execution/science/resource acceptance pending. The complete nine-path controller is frozen but unadopted and retains the parser defect; P30 must reconcile P01 and review the entire final patch. The flux-limited transport design is frozen and unapproved; P40 must settle the physical/domain decisions before source edits.

SpectralOperators 0.1.9 registration/tag/release are complete. HPIC requires that version; compatibility with actual minimum/latest Julia remains part of current/final qualification.

## File ownership and parallel work

1. Assign one integration lead for module/exports, dependencies, workflows, shared IO codecs/unions, full-test registration and global docs. Each actual path has exactly one writer across all worktrees.
2. Assign independent sizeable slices under the plan's parallel-lane table. Freeze equations, interfaces, aliases, publication, phase and persistence first. Propose shared-owner hunks instead of editing another writer's file.
3. Collision/source and transport paths overlap; controller/calendar and parser/current owners overlap; implicit/explicit EM and boundary operators overlap. Serialize those owners or transfer a frozen complete result explicitly.
4. Require a full patch/file handoff with source base/hash, complete caller review, schema, source checks and explicit pending gates. The lead reads the complete final result, reconciles current source and preserves canonical physical policies.
5. Run compact independent essential controls and actual full verification on GitHub/approved remote hardware. While waiting, advance another ready nonoverlapping slice and keep the tracker current. Source-only review never closes executable/science/resource gates.

## Shared-machine and publication rules

No local Julia/parser/compiler/formatter/test/numerical-proxy/profiler/MPI/device workload. Preserve unrelated dirty work and inherited processes/CI. Stop only an exact active-task process proved launched in the current session. Do not reset, stash or delete the original primary.

Keep all numerical assertions, compiler caps/filters, the 34-entry compiler ledger and fixed resource limits unchanged. Download and inspect raw failures; do not suppress a gate or classify actual code failures as billing. Actual hardware is required for backend claims.

Use the existing remote `Implementation essentials` workflow for focused suites. `suite=format-candidate` prepares an exact formatter artifact after the loader fix; review/apply its hunks only after the affected writer has handed off. Preparation success does not accept input formatting or physics. Re-run real full-CI formatting and complete exact-head checks afterward.

Push every coherent reviewed change to a validation branch. Promote that exact revision to main after all required CI passes; re-check any merge/rebase result whose source differs. Record immutable SHA, actual runtime/dependencies, raw logs and separate scientific/resource acceptance. Read back main and its CI after publication. Billing exceptions require direct evidence and never establish unexecuted science/hardware acceptance.

## Stop boundary and completion

The handoff plan covers 44 executable work packages and every one of the 79 original milestones. The 23 legacy-checked / 56-open count is 29.1% of an unweighted checklist; it does not measure full simulator completeness. All eleven capability areas retain work.

OpenCode resumes with P00/P01, then executes independent ready packages and dependency joins. Full closure requires the complete model/support matrix, all applicable G1–G10 evidence, final whole-code audit and measured optimization, current docs, minimum/latest Julia compatibility, accepted main publication and its green CI. Codex performs no further implementation after this handoff.
