# Distributed metadata-alias admission

## Pre-Edit Brief (2026-10-27, candidate 3)

- Mode: invalid-input robustness repair, with a portable test-only storage fixture.
- Requirement: enforce the existing read-only/private-workspace contract for detectable aliases of all twelve kfull/kvec/ik/midx metadata vectors before mutation, collectively on every rank. Preserve input/output APIs, arithmetic and communication.
- Baseline: HPIC765b233/SOd704f00. Preserved original two-rank kvec and candidate1 ik witnesses accept aliases and change metadata. Candidate2 adds the canonical eleven-line guard and passes initial both-runtime gates. Its direct pointer fixture assumes metadata bytes cover the full physical field, introducing a known Int32/Float64 size precondition. Julia1.13.3 officially supplies32-bit binaries; no32-bit runtime is available here. Preserve both earlier freezes and processes.
- Owner and invariants: `_require_physical_field` delegates every public operator's admission to `_aliases_workspace`. The grid's metadata are read-only and operands must not overlap private storage (`ext/pencil_fourier.jl`, `docs/DISTRIBUTED_DERIVATIVE_PLAN.md`). Rejection remains origin ArgumentError/peer ErrorException, before metadata/public/private-array mutation, followed by collective liveness. Undetectable unsafe aliases remain outside the supported adapter contract.
- Canonical strategy: preserve candidate2 production byte-for-byte: homogeneous real, complex and integer metadata loops in `_aliases_workspace`. No public method, formula, buffer, dependency or collective change.
- Portable fixture: retain ordinary reshaped real metadata. For different element types, allocate initialized byte storage large enough for both a logical metadata vector and complete physical field, including alignment slack. Create bounded own=false views at the same aligned start, copy original metadata values, and reconstruct the same concrete grid with that one numerically identical vector. The physical view may extend beyond the logical metadata span, but remains inside their shared owned backing allocation. Explicit GC.@preserve covers backing storage and metadata through all operations; no grid storage is freed. This removes the architecture condition rather than skipping Float64/midx. Compare metadata values and a valid derivative to the original grid before rejection tests.
- Compatibility: ordinary valid calls remain unchanged. All previously admitted invalid metadata aliases reject through the existing collective validation owner. Tests use only the same grid/field types, geometry, plan and private arrays; their custom backing allocation changes fixture ownership only.
- Verification: rerun the changed alias suite at actual MPI1/2/4/8 on Julia1.12.7/1.13.3; preserve candidate2 full operator/convolution, compiler/Aqua and same-dependency parent/candidate physical-byte/allocation controls because production, numerical tests and Project remain identical. Check fixture byte bounds/alignment, same values/type/pencil/dataids, all four families, fourteen operand surfaces, every origin, exact alias-error text, unchanged metadata/public/private sentinels and subsequent liveness. AST/idempotent formatter margin100, complete diff and source-pins review. Root owns integration/push and exact-SHA CI.
- Wrong patches rejected: identity-only or real-only checks; output-only or operator-specific guards; late scratch mutation; local-only throw; finite-input rejection substituted for alias admission; insufficient/misaligned backing storage. The measured aliases, all-family/all-origin matrix, private sentinels, error text and explicit span assertions reject these.

## Progress

- [x] Rebuild the contract and reproduce real and complex metadata admission failures.
- [x] Preserve canonical admission and make the fixture portable without an architecture skip.
- [x] Complete minimum/latest Julia MPI alias gates and unchanged-production compiler/resource controls.
- [x] Root registry-only MPI1/2 and compiler/Aqua integration on both runtimes.
- [x] Separate main push as `494faab`.
- [x] [All six exact-SHA CI jobs](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/40398004888).

All 38 required local commands terminated successfully on Julia 1.12.7 and 1.13.3: eight full operator/convolution jobs, four compiler/Aqua jobs, sixteen ordinary resource jobs, two paired resource follow-ups, and eight final portable-alias jobs. Both numerical and alias suites exercised actual MPI processes at 1, 2, 4 and 8 ranks. The final eight-rank alias suite passed 11,016 assertions per rank, plus 80 portable-span assertions per rank. The fixture explicitly checks Int32/Float64 storage bounds; actual 32-bit Julia execution remains unverified.

Candidate 3 production, other numerical tests and Project.toml are byte-identical to frozen candidate 2. Its compiler comparisons preserve all frames, signatures and messages in 22 complete report blocks per runtime after a checked source-line map; six reviewer mutations reject. Fresh registered SpectralOperators 0.1.5 resolution and the installed package tree were independently checked against tree `12831f35e5026cded6879554f3001dc1b614d01b`.

Ordinary controls cover 480 profiles and 1,440 allocation sample pairs per runtime. All original public physical-state hashes match the same-dependency parent, and no original candidate allocation sample increased. One original Julia 1.13.3 single-rank parent sample was higher than its candidate counterpart (22,576 versus 16,192 bytes); its cause remains unverified. One predeclared full 32-profile parent/candidate repeat matched all physical states and all three allocation samples exactly. Original measurements and the initial failed equality review remain preserved. These checks establish neither a timing improvement nor multi-node scaling.

The sealed evidence directory is `HybridPlasmaPIC-audit-evidence/2026-10-27/distributed-metadata-alias-v3/`, beside the package checkout. `final-local-audit.json` records the gate totals and limitations; `HANDOFF.md` identifies the integration commit and supporting receipts. Detection remains limited to aliases reported by `Base.mightalias`; this repair does not promise detection of hidden unsafe pointer overlaps.

Root integration adds eight terminal-zero commands on Julia 1.12.7 and 1.13.3.
Every MPI child verifies its active project, HPIC root and registered 0.1.5 tree.
The single-rank alias gate passes 2,244 assertions plus 80 portable-span checks;
the two-rank gate passes 2,856 plus 80 per rank. Twenty-four concrete analyses
and Aqua pass per runtime. Complete 29/27-report compiler aggregates and all
252/203 virtual-frame records match the frozen parent after exact root mapping.
The three executable changed files are byte-identical to sealed candidate 3;
rebasing retains the exact landed RK4 source. Root receipts are retained in
`takeover/metadata-root-integration`.
