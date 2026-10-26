# Distributed storage and derivative: Stage B pre-edit brief

- Mode: implement the bounded Stage B contract in DISTRIBUTED_OPERATOR_IMPLEMENTATION_PLAN.md.
- Requirement: collective Float32/Float64 real fields and first derivatives on owned periodic 3D pencils, with no production gather. Preserve serial mode arithmetic, global scaling, collective recoverable failure and delayed output publication.
- VERIFIED gap: current Pencil extension contains only DistributedPencilFFTPlan and transform/diagnostic methods; src/extensions.jl exports no distributed grid/field factory. SpectralOperators shared-kernel candidate owns the required metadata, scale, copy/store and derivative multiplier; its publication gate remains open.
- Canonical owners: factor the existing serialized configuration agreement inside the Pencil extension and use it from both factories. Existing distributed_fft_plan still owns all FFT construction. Shared SpectralOperators helpers own all Fourier arithmetic. The new HPIC grid owns its plan and three local work arrays.
- Contract: all ranks agree positive 3D integer shape and one common Float32/64 length type and value before topology. Checked global volume; exact physical pencil identity and Array-backed fields; common precision and logical axes. Each workspace belongs to one caller. Metadata vectors are read-only. Derivative axis is validated/agreed collectively; input/output may alias, but neither may overlap workspace storage. No communicator ownership transfer or field gather.
- Failures: local validation, metadata, scaling and multiplication errors use the existing collective stage; input nonfinite and lost-nonzero scaling reject before transposes. Stage complete real output in private scratch and publish only after global success. Fatal allocation, MPI process/transport failure remain existing exclusions. Unrepresentable mathematical outputs retain the serial floating-point behavior.
- Compatibility: serial API and FFT defaults/transform contracts unchanged. Add two extension factories and methods of the existing deriv!/deriv functions for the new HPIC-owned grid. Consumer compat will require the published shared-kernel dependency revision; local develop is verification only.
- Implementation: validate/agreed geometry; shared metadata; canonical FFT constructor; one complex physical, one complex spectral, one real scratch. Use the same physical backing order for real/complex copies and the actual global logical spectral view for multiplier indices. Precompute derivative gain from the shared metadata. Validate all operands in a local collective stage, agree axis, reduce input maximum, prepare scratch, transform/multiply/inverse/rescale, agree publication, copy to output.
- Verification: both runtimes, 1/2/4/8 ranks, odd/even unequal and thin shapes, all axes/precisions, analytic and independent dense DFT references; shared serial comparison; DC/Nyquist; huge constants; mismatched geometry/precision; wrong field ownership/type/axis and rank-local nonfinite/scaling-loss/workspace alias; unchanged failed outputs and later collective liveness. Full FFT regression, concrete/compiler-tree review, allocations/retained storage, explicit reductions/transposes and reproducible dependency resolution before promotion.
- Wrong patches: local indices for global modes; local volume/max scaling; missing inverse normalization; publishing before remote failure; rank-local dispatch errors; ignoring axis/length mismatch. Independent modes, dense DFT, uneven rank/layout and injected-origin rollback tests distinguish these. No timing claim without a separately frozen controlled campaign.

Implementation and results are recorded separately. This slice does not provide halos, deposition, particle ownership, a distributed stepper or checkpoint support.

## Candidate state

V1's standalone test omitted PencilFFTs and stopped at the existing missing-extension error before running the new code. The V2 test imports both declared extension dependencies. Its unchanged production passes the complete initial matrix at 1/2 ranks on Julia 1.12.7 and 1.13.2. V3 keeps those production bytes and integrates the test into the real CI entry, adding geometry/precision disagreement, invalid types, checked global-volume overflow, minimal valid pencils, extreme domain lengths, private-input aliases and distinct public wrappers sharing storage.

Landing (2026-10-26, on `1b7b8b7`): two corrections. `PencilArray` is an
immutable wrapper, so two constructions over the same storage compare `===`; the
alias test now asserts shared parent storage and `mightalias` instead of
`!==`. JET found three new reports because the collective-stage results were
untyped where `deriv!` is analyzed with abstract arguments; declaring
`j::Int` and `local_magnitude::Union{Float32,Float64}` (the precisions
`distributed_fft_plan` already enforces) restores the parent report set.
`test_distributed_fft.jl`, including this suite, passes at 1, 2, 4 and 8 MPI
ranks on Julia 1.12.7 and 1.13.2 with registered SpectralOperators 0.1.3;
`test_operators.jl`, `test_spectral_extra.jl`, Aqua and JET (24 / 21) pass.
Exact-commit CI is pending. Timing has not been measured.
