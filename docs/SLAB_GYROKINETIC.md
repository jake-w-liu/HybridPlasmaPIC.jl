# Full-f electrostatic slab gyrokinetics

This family evolves a self-consistent gyrocentre population and potential in a
3D periodic Cartesian box. Source implementation is separate from runtime,
scientific convergence, independent-code and resource acceptance. Those gates
remain open; the broader CAP-006 fusion/gyrokinetic/AMR goal is incomplete.

`GyrocentreSet{T}` stores three `X` vectors, `vpar`, physical `mu`, number
`weight`, `id`, `tag`, and fixed positive species `q,m`. It is a real full-f
population, with no artificial perpendicular velocity components. `T` is Float32
or Float64; `FourierGrid{3,T}` needs at least three nodes per axis. Both odd and
even grids are admitted. Each population may be empty or have zero-weight
tracers; the whole collection must carry positive full-f number.

```julia
g = FourierGrid((9, 9, 9), (8.0, 8.0, 8.0))
population = GyrocentreSet{Float64}(2; q=1, m=1)
population.X[1] .= (2.0, 6.0)
population.X[2] .= (3.0, 5.0)
population.X[3] .= (2.0, 6.0)
population.vpar .= (-0.1, 0.1)
population.mu .= 0.01
species = [population]
simulation = SlabGyrokineticPIC(g, species; B0=1, Te=1, e=1, nring=16)
init!(simulation, species)
step!(simulation, species, 0.001)
report = slab_gyrokinetic_diagnostics(simulation)
save_checkpoint("slab.bin", simulation, species)
load_checkpoint!(simulation, species, "slab.bin")
```

The example is an API illustration, not a validated physical equilibrium or
accepted numerical result. Choose the timestep and marker sampling through the
physical problem and convergence studies; there is no gyrokinetic timestep
certificate borrowed from the hybrid whistler recommendation.

## Physical model and units

The fixed field is `(0,0,B0)`, with B0 positive. The model is collisionless,
electrostatic and has uniform equilibrium densities. It has no mirror,
curvature, toroidal geometry, nonlinear polarization, kinetic electrons,
magnetic perturbation, heating/source, transport, MPI, device or adaptive mesh
integration. Magnetic moment and full-f weights are fixed; changing those
invariants, species q/m or particle counts requires a new object.

Equations use coherent SI form. Te is an **energy**, not an electron-volt count
or temperature in kelvin. The default e=1 is appropriate only when the selected
charge unit gives that value; supply a consistent elementary-charge magnitude
for SI data. There is no implicit Gaussian c, vacuum permittivity or hybrid
mu0 normalization. For normalized inputs the reference scales must satisfy
`v*=L*/t*`, `B*=m*/(q*t*)`, `phi*=m*L*²/(q*t*²)`, and
`mu*=m*L*²/(t*²*B*)`, with consistent number/volume units. Archive metadata defaults
to this coherent SI-form convention; conversion through PlasmaUnits is not
performed.

Stored mu means `m*vperp²/(2B0)`. Thus
`rho²=2*m*mu/(q²*B0)`. Radius construction retains the complete source factors
through the shared exact square-root owner. Positive mu must yield a finite
nonzero radius. Zero mu has exactly zero radius, with the same ring count.

Physical number weights already include their initial phase-space quadrature
and Jacobian. Let `W_s=sum(weight_s)`, `V=prod(L)`, `n0_s=W_s/V`. The exact
background is `rho0_Q=sum(q_s*W_s)/V`, `ne0_Q=rho0_Q/e`. It neutralizes the exact
inventory. Rounded ne0_display is metadata, never the owner of that identity.
The field uses independently rounded susceptibilities
`chi=sum(m_s*n0_s)/B0²`, `alpha=e²*ne0_Q/Te`; both and displayed ne0 must be
positive finite normal values.

Define P as subtraction of the y-z plane average at each x. The electron law is
`ne=ne0_Q+(alpha/e)*P(phi)`. It has no response to zonal y/z-zero modes. Every
stage requires its complete electron density to materialize finite and positive;
there is no density floor. This is a linear modified adiabatic approximation.
The reported response amplitude and k-perp-radius coverage do not certify the
small-fluctuation/long-wavelength ordering.

The equilibrium polarization and modified adiabatic ingredients follow the
model distinctions in [Lanti et al., ORB5, sections 2.2, 2.3.1 and 2.3.3](https://arxiv.org/pdf/1905.01906v1).
Their toroidal discretization and validation do not transfer to this slab family.
The physical magnetic-moment convention is also given by
[Fitzpatrick, guiding-centre motion](https://farside.ph.utexas.edu/teaching/plasma/Plasmahtml/node18.html).

## Coupling, field equation and time phase

Every fixed-radius ring uses the existing robust geometry and phase order. The
existing TSC basis deposits physical gyroaveraged ion charge and interpolates
the same potential. Force differentiates that interpolated basis, using centered
potential differences that annihilate a constant represented potential. It is
not a TSC gather of a separately differentiated spectral electric field.

The field equation is `chi*(-Delta_perp)phi+alpha*P(phi)=Pi0(rho_gc)`, with zero
potential DC. Its symbol uses full perpendicular Fourier k, including even-grid
Nyquist. Only DC is null; using the centered first-derivative k squared would
create an incorrect unscreened zonal checkerboard nullspace. The non-DC symbols
must be finite positive normal values with finite inverse gain.

A fresh potential transform checks the original coefficient expression against
the source transform before acceptance. The norm is the maximum real/imaginary
modal residual divided by grid count, with `atol+rtol*source_norm` in density
units. Gauge allowance is `atol/lambda_min+rtol*potential_norm` in potential
units. These comparisons retain exact physical scales; saturated displayed
receipts cannot cause Inf/Inf acceptance. Source DC mean, charge-deposition
defect and gauge residual remain visible.

Characteristics are `Xdot=(-psi_y/B0,psi_x/B0,vpar)` and
`vpardot=-(q/m)*psi_z`, retaining original charge/mass factors. All four classical
RK4 stages re-deposit the entire collection and solve the common field. A fifth
endpoint refresh synchronizes stored phi/rho with completed X/vpar. The existing
RK affine owner retains original durations through range repair. No field is
frozen across stages, and no velocity is at a hidden half time.

Initialization privately wraps initial coordinates and constructs the field,
then publishes time=step=last_dt=0. Actual zero dt is a validated no-op;
nonzero inputs disappearing under conversion reject. Positive steps require a
representably advancing clock. Failed stages leave all public physical arrays,
clock and previous receipt unchanged. There is no automatic timestep retry or
partial species publication. Public arrays/owners retain identity and count.

Nearest rounding and gradual underflow are checked entry conditions. Lost
nonzero FFT rescaling or modal division rejects before publication; the family
does not promise arbitrary active-mode FFT cancellation accuracy. Ordinary
CPU storage and all public/private/FFT/metadata roles must be independent,
including empty arrays. Give each concurrent simulation its own grid and work.
The busy guard prevents reentry; it is not synchronization against external
mutation. Ready public arrays are read-only between calls; use init! after
editing positions/parallel velocities or fields.

## Diagnostics and restart

Diagnostics report exact number/charge inventories, parallel momentum,
parallel kinetic energy, fixed mu*B energy, equilibrium polarization energy,
quadratic electron-response free energy and their complete total. The field
terms are computed from the full-k spectral quadratic form with the correct
FFT normalization. A signed source-pairing defect exposes numerical disagreement.
Saturated scalar outputs carry range flags. No empirical marker-log-weight
entropy estimator or exact finite-step energy/momentum theorem is claimed.

Checkpoint schema 1 `:slab_gyrokinetic_fullf` saves every population leaf, phi,
rho_gc, synchronous clock, preceding interval, scalar field receipt, actual
quadrature table and exact fixed configuration policies. A destination must
match species order/count/q/m/mu/weights, geometry, precision and field controls.
Original incoming aliases are checked before whole-snapshot detachment and
revalidation. Incoming density is privately re-deposited and saved potential
passes the same fresh field/gauge/electron gates. These validation computations
have a cost; the saved fields still restore verbatim, without initialization,
a replacement solve or recentering. Stale destination positions/velocities,
fields and clock can be replaced. Full fixed-length physical publication occurs
only after all validation and preparation.

`archive_run` uses the same physical payload and shared checksum container.
Its required rng_seed argument is provenance metadata; this deterministic family
has no active RNG state. `load_archive!(simulation, species, "run.bin")` first
verifies the shared run-container schema/checksum, then uses the same complete
physical restore as `load_checkpoint!`; omitting species selects the bound owners.
`load_archive(path)` remains a read-only data/provenance reader. A malformed
container or rejected destination leaves the public physical phase unchanged.

## Acceptance status

`test_slab_gyrokinetic.jl` covers boundary wrapping, the zonal full-k
polarization solve against an analytic amplitude, gather consistency, the
zero/nonpositive-dt guards, the parallel free-stream characteristic, uniform
equilibrium deposition, checkpoint and archive round-trips and the diagnostic
inventory on Float32/Float64.

`test_slab_gyrokinetic_oracles.jl` adds physics oracles on Julia 1.12.7 and
1.13.3, Float32 and Float64: an independently coded ring-TSC deposit (three
markers including finite and zero `mu`) matching `rho_gc` cell by cell and the
stored radii matching `sqrt(2 m mu/(q^2 B0))`; an independent FFTW solve of
`chi*(-Delta_perp)phi + alpha*P(phi) = rho` agreeing with `st.phi`, the stored
symbol table `st._lambda` and a recomputed modal residual inside the receipt
tolerance; per-marker E×B drift `-psi_y/B0` and parallel `vpar` transport on a
y-structured sheared configuration, where uniform x/z lattices keep `rho`
invariant so one RK4 step is an exact translation checked marker by marker;
and independent quadratic-form recomputation of the polarization, electron
free-energy, pairing-defect, parallel and magnetic-moment reservoirs. Guard
coverage includes weight/mu invariant mutation rejections with serialized
state intact and constructor domain failures (negative `mu`, `B0 <= 0`,
`Te <= 0`, `nring < 3`).

Still open for 006B acceptance: genuinely nonuniform-time (RK truncation)
order measurement, an analytic linear-mode dispersion oracle, mirror/bounce
orbits (absent in this constant-B model by contract), delta-f/weight-evolution
variants, field-solve failure rollback through a stage failure rather than the
boundary guards, independent-code (e.g. GX-class) comparison, marker/ring/field
resource budgets and MPI/device coverage. Runtime continuation, accuracy/
refinement, package/JET, latest Julia and cost acceptance require
exact-revision remote evidence.
