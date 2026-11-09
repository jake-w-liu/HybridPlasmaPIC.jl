# Frozen OpenCode handoff packet — 2026-11-09

Start with the [implementation plan](../../docs/OPENCODE_IMPLEMENTATION_HANDOFF.md) and [79-row acceptance tracker](../../docs/OPENCODE_ACCEPTANCE_TRACKER.md). Codex stops after preparing this packet. OpenCode has not started and no Codex writer remains active.

- [MANIFEST.json](MANIFEST.json): exact code/branch/main snapshot, protected primary hashes, historical compiler ledger and SHA-256 inventory. Manifest itself is protected by the containing Git commit.
- [CODE_TREE_BASELINE.txt](CODE_TREE_BASELINE.txt): Git tree listing for `src`, `ext`, `test`, `Project.toml` and `.github/workflows` at `5aa628513fb854c96b2b387780edef1ea1faa66b`. Final handoff source must match it.
- [WORK_PACKAGES.json](WORK_PACKAGES.json): 44 packages, real owner paths, acceptance dependency DAG, deliverables and closure rules.
- [ROADMAP_CROSSWALK.json](ROADMAP_CROSSWALK.json): all 79 immutable milestone texts, line references, hashes and package mappings; 23 legacy checked, 56 open.
- [ISSUE_REGISTER.json](ISSUE_REGISTER.json): verified loader/JET/format/focused failures and explicitly unexecuted hypotheses. First action is P00/P01; the current source is not package-load accepted.
- [controller final receipt](controller/FINAL_RECEIPT.json) and [complete patch](controller/coupled-empic-final-v1.patch): nine frozen paths based on `566ddc1`, mechanically applicable but unadopted, final source review and execution pending. It retains the two parser errors; reconcile P01 before P30 adoption. No standalone Julia source from this packet enters compilation.
- [flux proposal](flux-proposal/PRE_EDIT_CONTRACT.md), [primary sources](flux-proposal/PRIMARY_SOURCES.md) and [remaining source/transport slices](flux-proposal/REMAINING_CAP010_CAP011.md): frozen read-only design, no approval or implementation. Preserve the explicit pending design choices.
- [CAP-002–006 review](reviews/CAP002_006_REVIEW.md) and [CAP-001/007–011 review](reviews/CAP001_007_011_REVIEW.md): finite read-only planning audits. Their source snapshot was `5de6522` plus then-uncommitted CGL work; final CGL source is committed at `5aa6285`. The main plan's exact starting-state table supersedes historical assignment/state prose.
- [CGL source checks](cgl-reviewed/ROOT_SOURCE_CHECKS.json), [root review](cgl-reviewed/ROOT_SOURCE_REVIEW.md) and [refinements](cgl-reviewed/ROOT_REFINEMENTS.json): checked source/text scope only. The worker-original patch is provenance, not the final integrated patch.
- [CI snapshot](ci/STATUS_BEFORE_HANDOFF_PUSH.json), [new failure summaries](ci/LATEST_FAILURE_LOG_SUMMARIES.json) and [MPI artifact receipt](ci/MPI_FAILURE_ARTIFACT.json): inspected historical exact-head observations. New handoff-branch jobs are not accepted by these earlier results. Retained loader/MPI raw logs accompany their receipts where available; full compiler reports can be re-downloaded from the named GitHub jobs.

Run from the repository root, using only standard-library Python and read-only Git:

```sh
python3 handoff/2026-11-09/verify_handoff.py
```

It checks every packet hash, the active handoff documents, complete roadmap mapping, dependency DAG, actual owner paths and the original source fingerprint. This is a metadata check, not a Julia parser/runtime/science/performance test. After source-only edits, use `--allow-source-changes` to check the immutable packet without requiring the new source to equal the historical baseline. Active docs and trackers must evolve during implementation; once those change, verify this initial handoff in a clean checkout of its published commit. Do not rewrite historical manifest hashes to hide changed evidence.

Frozen reports retain original absolute paths, sibling worktree references and historical line numbers as provenance. They are immutable historical records; use the indexed local copies and current main plan for operation. Do not treat an old writer assignment, hypothesis, approval of one bounded contract or an old green CI run as current implementation acceptance. Preserve hashes when quoting evidence; add a dated superseding record instead of rewriting historical receipts.
