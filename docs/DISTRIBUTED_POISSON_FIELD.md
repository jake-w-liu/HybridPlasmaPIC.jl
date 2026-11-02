# Rank-owned periodic Poisson field

`distributed_poisson_field!(E, rho, grid; neutrality=:require)` supplies a static
3D CPU Float32/64 electrostatic field solve on the existing physical pencils.
It reuses the canonical serial Poisson mode map and the existing distributed
transforms. The [distributed ES collection](DISTRIBUTED_ELECTROSTATIC_IMPLEMENTATION.md)
composes it with actual-pencil ownership, complete CIC transport, migration and
the shared leapfrog step. Broader distributed families and runtime acceptance
remain open.

Import PencilArrays and PencilFFTs with HybridPlasmaPIC, construct a
`distributed_fourier_grid`, and allocate rho and three independent E components
with `distributed_field`. Fill owned cells through `PencilArrays.global_view`.
Call the solver on every grid rank in the same order. Grid n/L/stored spacing
and neutrality policy must agree; metadata remain read-only. Nearest rounding
and gradual underflow are required. No tests or runtime execution are added
during this implementation-first pass; acceptance is deferred.

## Source, null modes and charge ledger

The input is finite stored charge density, with field convention `div(E)=rho`
on derivative-resolved modes. The shared scalar map is `Ehat=-im*k*rhohat/k²`.
It sets every derivative-null DC/Nyquist combination to zero, matching the
existing serial electrostatic solver. No complete Gauss identity for unresolved
charge modes or arbitrary FFT roundoff is claimed.

The exact sum of all represented rho entries is reduced as scalar data, with
no premature rank-local IEEE rounding. `:require` accepts only exact zero total.
`:neutralize` permits the explicit uniform compensating charge; the Poisson
null-mode policy supplies no uniform electric field. The returned NamedTuple is

```text
(E = output_components,
 source_charge = exact_sum(rho)*prod(exact_stored_dx),
 neutralizing_charge = :neutralize ? -source_charge : 0)
```

Charges are exact Rational{BigInt} values. The original density is not altered,
unless it is explicitly also one output. This function does not infer particle
weights, background species, physical count conversions or boundary charge.

## Transform and publication contract

One output component may alias rho. Output components must be mutually
independent and use the grid's exact physical pencil, precision, axes and
ordinary CPU backing, disjoint from private Fourier storage and metadata.
All source data and modes are consumed before output publication.

The solve uses one forward and three normalized inverse transforms, followed
by the existing inverse rescale. The component multiplier is bounded by
`1/k_min` over nonzero derivative-resolved axis entries. This bound rounds
upward before selecting the existing global-count FFT scale. An unrepresentable
bound rejects conservatively. Lost-nonzero input rescaling and the existing
finite-coefficient/FFT limits remain; this is not an original-source range or
cancellation-relative accuracy completion. A completely derivative-null grid
fills private zero fields directly after charge admission.

Each local fallible stage agrees success before the next communication. All
three finite output candidates and the exact return ledger are prepared before
any public store. A recoverable rejection preserves output components; private
FFT scratch may change and is overwritten on reuse. The caller owns the
communicator and must exclude concurrent use of the grid or arrays. MPI
process/transport loss and fatal resource failures retain existing exclusions.

The scalar reduction is shared with distributed Maxwell energy. It exchanges
only exact scalar payloads in rank order; no density or field is gathered.
Coefficient/array reuse is implemented, but runtime, compiler/format, actual-rank
failure/liveness, numerical and memory/cost acceptance remain unverified.
