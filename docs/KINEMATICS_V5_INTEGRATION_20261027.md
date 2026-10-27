# V5 integration onto b583e41

## Pre-edit brief

- **Mode:** integrate the sealed Lorentz-gap/priming repair; mechanically move
  the landed radiation arithmetic to the canonical EFT utility.
- **Requirement:** retain V5 executable kinematics and all three companions;
  retain the complete radiation expression landed at
  `07792ab5b7ac5854948fe85fcb6ab196d43a3366`, changing only the mapped private
  arithmetic symbols and their definition locations. No V6 optimization.
- **Baseline evidence:** this new worktree starts clean at
  `b583e41f99e833cfc407ad764c1303d319b41c7e`. `git diff cd55831 b583e41`
  shows no changes in the six pre-existing V5 kinematics/module/test owners;
  radiation and `Project.toml` have independent landed changes. The radiation
  source and its two regression files are byte-identical between 07792ab and
  b583e41. The sealed V5 patch applies cleanly with radiation excluded.
- **Design gap:** replacing radiation with V5's old-base file would discard
  the exact-invariant completion and its amplification guard. Its pure EFT/DD
  primitives still have the original bodies that V5 extracted.
- **Contracts/invariants:** preserve current radiation arithmetic, admission,
  exact cold tier, task-local state, validation and normwise result contract;
  preserve V5's local kinematics contract and documented compatibility slices.
  The representable output cone and trajectory-wide contract remain separate.
- **Canonical owners:** shared primitives belong in
  `src/utils/error_free_arithmetic.jl`; radiation retains all physical policy
  in `src/particles/radiation_reaction.jl`; V5 retains its existing Boris and
  priming owners. No compatibility wrappers or second arithmetic owner.
- **Compatibility:** apply only the sealed private-symbol mapping, move each
  original complete definition once, preserve its executable AST, and preserve
  every other radiation expression. Current dependencies and root's other
  landed changes remain intact.
- **Strategy:** apply the checked V5 patch excluding radiation, then the three
  companions. Extract the original primitive definitions from current
  radiation and rename their remaining references with the sealed map. Prove
  complete radiation reconstruction modulo this bijection, with expression,
  definition-removal and arithmetic-mutation controls. Compare all other
  changed executable trees to frozen V5; boundary test additions are explicit.
- **Verification:** prepare source/harness/environment pins first. Coordinate
  with root before expensive tests. Then run the unchanged 54,877 affected
  radiation assertions, 108 cold resource/task checks and 11,447 independent
  checks on Julia 1.12.7/1.13.3, plus integrated quality and V5 checks needed
  for changed dependencies or source. Reuse frozen receipts only with actual
  source/AST equivalence, keeping their original dependency scope explicit.
  Never replace a failed gate with a weaker gate.

## Progress

- [x] Create clean isolated worktree at b583e41; preserve existing worktrees.
- [x] Read landed radiation policy, V5 patch/companions and existing harnesses.
- [x] Integrate source and prove complete expression equivalence (38 assertions
  on both runtimes after formatting; full radiation reconstruction plus mutation
  controls and unchanged V5 numerical owners).
- [x] Prepare six pinned registry-only SO 0.1.6 environments and serial commands;
  verify the installed Git tree and complete declared Manifest changes.
- [x] Pass fresh radiation gates on Julia 1.12.7/1.13.3: 54,877 affected
  assertions, 108 cold resource/task checks and 11,447 independent checks each.
  All 2,432 raw independent cases match bitwise; all 21 cold allocation profiles
  match the original root repair exactly.
- [x] Complete fresh permanent kinematics, actual-step, task, resource and Metal
  gates on both runtimes: original oracle 10 + actual steps 8; Metal 824; caller
  contexts 512 each; 24 profiles/72 warm allocation samples all zero. The
  standalone permanent wrapper passes 39,658 assertions each, including an
  accidental repeat of the 18,009-assertion factor block; no extra coverage is
  claimed for that repeat. The canonical suite includes it once.
- [x] Complete current-dependency baseline/candidate Aqua and compiler gates:
  19 concrete calls clean per runtime and worktree; all 29/27 inherited package
  report trees unchanged after checked source/root/private-symbol mapping.
- [x] Review task-owned source, documentation, validation and complete diff;
  prepare the reviewable handoff and apply-check it against clean b583e41.
- [x] Root: merge into current main and push as `fd0a168`.
- [x] [All six exact-SHA GitHub CI jobs](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/40407584532).

Frozen V5 evidence and its limitations remain in
`KINEMATICS_AUDIT_20261027.md` and `PARTICLE_PRIMING_INTERVALS.md`.

The registry-only dependency is SpectralOperators 0.1.6, tree
`ea7046c76096424f364ea15a00c1a32a9bf3f6db`. A complete official General
snapshot at commit `d27db9cf6493fd601c486a4cc6a5b68825b47648` supplies the
entry in a task-local depot because the package-server snapshot initially
exposed only 0.1.5. No dependency path/repository override is used. Julia 1.13
records the two reachable General names as `[General, General]`; this exact
metadata migration is checked separately from versions and source trees.
The 741-file freeze is checked around every numerical/quality command.

## Root integration

The source is merged onto `36d1794` in an isolated root worktree. All eleven
changed source/test owners remain byte-identical to the sealed candidate.
The landed midpoint, distributed metadata and BigFloat/RK4 owners are retained.
The root campaign passes all 23 commands on Julia 1.12.7 and 1.13.3:

- 21,649 permanent assertions, with the factor block included once; 3,835
  midpoint assertions and the original near-light oracle/actual first-step gates.
- 54,877 affected radiation assertions and 108 cold resource/task assertions;
  all 21 cold profiles retain their sealed allocation counts.
- Twelve warmed kinematics profiles retain exact output bits and 0 B over
  three samples; 512 caller-context calls pass per runtime. No latency claim.
- 1,240 assertions on actual functional Metal, including midpoint integration.
- Aqua and 48 concrete compiler calls pass. Full 29/27-report forests retain
  every message and actual frame field after checked unchanged-source line
  movement and seven exact generated-label occurrences. Eighteen latest-runtime
  Serialization frame paths resolve to the same complete release source;
  Julia's own path resolver and whole-file hash prove that relocation.

The strict compiler comparison was held before that path proof, and its diff
is retained. Twenty-five changed-message/frame/path/signature controls reject.
The root evidence is in `takeover/kinematics-root-integration/` in the dated
external audit archive. `fd0a168` is pushed to main; exact-SHA CI remains.
