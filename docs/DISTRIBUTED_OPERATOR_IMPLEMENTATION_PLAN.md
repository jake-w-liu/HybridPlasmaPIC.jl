# Domain-owned periodic operators: CAP-003A/B

Design source review: 2026-10-26. HPIC `00f7922`; sibling SpectralOperators
`cc515cf`, version 0.1.2. This is an implementation contract, not an accepted
operator or solver. The main completeness plan remains the scope owner.

## Checked architecture and pre-edit brief

- Mode: implement distributed physical storage and the periodic operator family.
- Requirement: real Float32/Float64 fields on fully periodic 3D pencils; derivative,
  gradient, divergence, curl, Laplacian, solenoidal projection and existing filters
  use owned cells and planned transposes. No full-field gather in these operations.
  Serial FourierGrid remains the numerical and compatibility reference.
- VERIFIED source gap: HPIC's Pencil extension exposes transforms and a diagnostic,
  but no distributed field operator. `HybridFields` and existing stepper workspaces
  retain global Array storage; this stage does not imply a distributed stepper.
- Canonical owners: SpectralOperators owns signed modes, odd-derivative Nyquist
  handling, range scaling and mode arithmetic. HPIC's optional Pencil extension
  owns MPI layout, validation agreement and transform transport. Existing physics
  methods and serial storage remain unchanged by this stage.
- Compatibility: retain all current serial signatures, alias rules, operation
  order and range behavior. Add distributed methods to the same imported operator
  functions on a new HPIC-owned grid/workspace type. No type piracy, duplicate
  Fourier physics or mandatory Pencil dependency in either core package.
- Strategy: extract only the existing mode/metadata kernels actually required by
  both backends in SpectralOperators, with serial differential checks. Introduce
  a concrete distributed workspace and orchestration in the HPIC extension.
  A generic grid hierarchy is not needed for this bounded reuse.
- Verification: frozen serial source, independent analytic modes, actual logical
  pencil indices, 1/2/4/8 real ranks and both Julia versions. Include rank-local
  invalid inputs, asymmetrical partitions, finite extreme amplitudes, DC/Nyquist,
  alias behavior, workspace reuse and exact ordinary allocation comparisons.

Checked sources: `SpectralOperators/src/periodic_fourier.jl` constructs
`midx/kfull/kvec/ik` and owns `_fft_safe_magnitude`, `_fft_input_scale`,
`_copy_fft_input!` and output rescaling. `differential_operators.jl` owns the
first-derivative, curl and Laplacian formulas; `divergence_control.jl` owns the
normalized projection; `filters.jl` owns the mode-ratio filter and strict 2/3
cutoff. Source hashes accompany this design in the audit evidence directory.

`PencilArrays/src/global_view.jl` implements global *logical* indices through an
OffsetArray around a PencilArray and accounts for permutation. Its view does not
copy. `PencilFFTs/src/plans.jl` accepts a real precision type as the third argument
of the Pencil constructor. Its forward output may have a different distribution
and physical order from its input. `pencil_input` and `pencil_output` expose those
layouts. Multipliers must therefore consume the output's logical global indices.

## Public storage and lifecycle contract

Proposed entry: `distributed_fourier_grid(n, L; comm=MPI.COMM_WORLD)` with matching
Float32 or Float64 domain lengths, plus a constructor helper that allocates a
real field on the workspace's physical input layout. Final names are fixed before
production edits and recorded in the next pre-edit brief.

Use an out-of-place full-complex FFT first, matching the serial mode convention.
Keep one reusable physical complex input, reusable spectral component buffers and
inverse scratch, plus physical real output scratch needed for complete publication.
Allocate through PencilFFTs/PencilArrays using the actual plan pencils. Use global
logical views for multiplier loops; do not infer permuted indices from raw parent
array dimensions. Replicated per-axis metadata costs O(nx+ny+nz); physical and
spectral fields cost O(owned cells). Measure both, including plan transpose buffers.

One workspace is used by one task at a time. Its communicator is caller-owned;
construction does not finalize MPI or free the caller's communicator. No global
mutable cache. Plans and workspace arrays have the same lifetime. Inputs/outputs
must use the compatible physical pencil and precision; matching local dimensions
alone do not establish ownership compatibility.

The first storage stage has no particle halos. It records owned physical ranges,
spacing, volumes and global coordinates; shape-radius halo exchange belongs to
003C before any local deposition/gather claim. Unsupported empty pencils reject
collectively, as the accepted transform constructor already specifies.

## Canonical numerical contract

1. Extract shared Fourier metadata construction without changing serial formulas:
   signed integer mode, 2π/L scale, finite spacing/wavenumber checks, full Nyquist
   for Laplacian/filter and zero Nyquist for odd derivatives/projection.
2. Reuse the existing first-derivative/curl multiplier kernels on globally indexed
   views. Extract the existing pointwise Laplacian, normalized projection and
   filter kernels once and call them from both serial and distributed owners.
   Keep each operation's arithmetic order; no new independent formula in HPIC.
3. Before each transform family, obtain the maximum finite input magnitude across
   ranks and components. The safety bound uses the **global transform length**,
   not local volume. All ranks use the same scale. Reject any nonzero value lost
   to rescaling on any rank before the next transpose. Use the existing scalar
   scale rule and representability guards.
4. Forward is unnormalized; inverse uses the dependency's normalized `ldiv!`.
   Apply physical rescaling exactly once after inverse. Preserve the serial
   projection's derivative-resolved null modes, including mixed DC/Nyquist modes.
5. Validate complete output/component alias constraints before communication.
   Scalar derivative/Laplacian may alias their input; curl outputs remain distinct
   from inputs and each other. Projection components remain pairwise distinct.
   Snapshot only where the existing valid alias contract requires it.
6. Stage output publication after successful shared validation and numerical work.
   Recoverable rank-local validation failures must give every rank an error before
   another collective. Unsupported MPI transport/process failures remain explicit
   exclusions; do not turn a transport error into an apparently successful result.

A bounded constructor prerequisite adds precision selection to the existing FFT
plan API, includes it in rank agreement, and passes it through to GlobalFFTParams
and PencilFFTPlan. Default Float64 behavior is preserved. Float16/BigFloat and
rank-mismatched precisions reject explicitly. This is a separate tested change.

## Stages and fixed gates

| Stage | Implementation | Acceptance evidence |
| --- | --- | --- |
| P | Collective FFT precision | Default compatibility; Float32/64 transforms and storage; rank mismatch and unsupported types; all family/rank matrices |
| A | Shared metadata and pointwise kernels | Serial complete-state differential, all public regressions, exact allocation comparison, compiler and both-version CI; publish a dependency revision consumed reproducibly by HPIC |
| B | Owned real storage and derivative | Uneven shape/rank ownership, global-coordinate modes, DC/Nyquist, in-place output alias, finite huge constants, one-rank invalid values, independent serial/dense-DFT references |
| C | Gradient/divergence/curl/Laplacian | Same-source kernels; divergence of curl, integration by parts, mixed Nyquist, Laplacian underflow/overflow behavior; output alias and failure publication |
| D | Projection and filters | Idempotence, divergence removal, nullspace retention, scaled k/field projection, strict 2/3 boundary, singletons and mode-ratio filtering |
| E | Resource and integration closure | Reused plans, bounded allocations/storage, full reports, transport inventory, docs/examples and existing six-job CI; no hot-path global gather |

For each stage, freeze source and the full workload matrix before execution. Cover
odd/even unequal axes, all physical directions, process grids 1/2/4/8 and both
precisions on minimum/latest Julia. Use serial equality only where arithmetic and
FFT ordering allow it; scientific error/order bounds are specified before the
run. Independent dense DFT and analytic references prevent circular shared-kernel
validation. Reject wrong local-index, wrong Nyquist, local-only scaling, omitted
inverse normalization, alias overwrite and dropped-rank-error mutants.

Timing is a separate frozen workload gate with parent/control measurements in an
uncontended window. Count global reductions, transposes and maximum rank-local
storage explicitly. A single-machine MPI result does not establish multi-node
scaling or GPU support. CAP-003C/D still require halo moments, particle ownership,
real hybrid/ES/EM evolution, collision composition, checkpoint and restart tests.

## Progress

- [x] Inspect current owners, public signatures and dependency layouts.
- [x] Specify shared numerical ownership and staged integration contract.
- [x] Pin/review precision prerequisite and serial kernel extraction briefs.
- [x] Stage P (precision) landed as `3a9c85d`; stage A (shared kernels) released
  as SpectralOperators 0.1.3; stage B (owned storage and derivative) landed
  after `1b7b8b7`.
- [ ] Implement and verify stages C/D/E separately; push each accepted slice.
- [ ] Integrate local field ownership with 003C and existing solver workspaces.

## Stage B concrete implementation boundary

The first operator slice fixes these public names before code changes:
`distributed_fourier_grid(n, L; comm=MPI.COMM_WORLD)`,
`distributed_field(grid)` and the existing `deriv!`/`deriv` function owners.
The grid factory and field allocator use the existing optional-extension hook
pattern. Derivative dispatch adds the HPIC-owned concrete grid type to the
imported SpectralOperators function; it does not redefine a dependency type.

Grid lengths must be a three-tuple of one common Float32 or Float64 type, with
positive finite values and representable spacing/wavenumbers. Shape conversion
uses the existing checked integer owner. Precision follows the lengths. Agree
all geometry and precision collectively before constructing topology; factor
the existing configuration-agreement transport into one private owner if it
is needed by both factories. Do not copy its serialization/error protocol.
The global transform element count is checked before using it in scaling.

Initially own one physical complex pencil, one spectral complex pencil and one
physical real publication scratch array. Derivatives can multiply the spectrum
in place through the shared safe exact-object alias contract. Stage C adds only
the extra component storage its complete vector operations actually require.
Each field allocator returns a zeroed CPU PencilArray using the exact physical
pencil object owned by the grid. Ordinary Array backing is the first supported
storage; GPU residency and arbitrary array adapters are separate stages.

The derivative validates all arguments collectively, including integer axis,
precision, exact pencil ownership, local axes and disjointness from private
workspace storage. Input and output may alias each other. `Base.dataids` is
implemented by PencilArrays, so storage aliases must be checked through the
actual arrays rather than object identity alone. All ranks must agree the axis.
A mismatched but individually valid axis must fail before the forward transform.

Compute the input magnitude locally during validation, reduce its maximum across
ranks, and call the shared scale rule with the global element count and global
wavenumber gain. Prepare only workspace storage; agree any lost-nonzero scaling
failure before the first transpose. Apply the shared multiplier to global logical
views, agree local multiplier errors before inverse transport, and rescale into
the real publication scratch. Publish public output only after the last required
collective success check. MPI transport/process loss and fatal memory errors keep
the existing explicit exclusions. No new guarantee that an unrepresentable
mathematical output is finite is inferred from the serial implementation.

The bounded gate matrix includes odd/even/asymmetric/thin grids, both precisions,
all axes and 1/2/4/8 ranks; analytic modes plus an independent dense DFT on small
grids; DC and mixed Nyquist; huge finite constants; global-versus-local scaling;
input/output aliasing; one-rank wrong ownership/type/axis/nonfinite values and
workspace aliases; subsequent collective liveness and unchanged rejected output.
Reuse the same grid repeatedly and record reductions, transposes, local storage,
ordinary allocation and complete compiler trees. No field gather is allowed in
production; test-only gathered serial comparisons are explicitly outside the
operator. This stage supplies no particle halo, current deposition, coupled
stepper or checkpoint compatibility claim.
