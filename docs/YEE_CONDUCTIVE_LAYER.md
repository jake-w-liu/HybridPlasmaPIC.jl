# Graded scalar-conductivity Yee absorber

Status: source candidate; no tests were added or run. Reflection, numerical and
long-time behavior, actual rollback/restart, compiler/format and complete-step
resource acceptance remain unverified. This is a field-only part of CAP-002C,
with no claim that bounded PIC or radiating-boundary work is complete.

`YeeConductiveLayer` surrounds the existing synchronized `MaxwellVacuum` KDK map
with matched scalar electric/magnetic conductivity evolution. It supports serial
CPU Float32/64 in D=1–3 and all three vector components. This is an artificial
graded absorber, **not PML/CPML or an exact outgoing boundary**. The
[Meep Absorber documentation](https://meep.readthedocs.io/en/latest/Python_User_Interface/#absorber)
makes the same scalar-conductivity/PML distinction. The separate
[`YeePML` source](YEE_PML.md) supplies causal directional auxiliary evolution and
persistence; its reflection/stability and runtime acceptance remain unverified.

```julia
vacuum = MaxwellVacuum(g; c)
# Seed vacuum.fields at synchronized time; explicitly impose valid PEC data.
layer = YeeConductiveLayer(vacuum;
    dt=h, cells=((nlow, nhigh),),
    peak_rate=((rate_low, rate_high),), profile=:quadratic)
init!(layer)
step_maxwell!(layer, h)
d = absorber_diagnostics(layer)
save_checkpoint(path, layer)
archive_run(archive_path, layer; rng_seed)
```

The pair tuples above illustrate D1; supply one pair per axis in D2/3. The
supplied vacuum must be unadvanced. Its physical fields/clock and candidate are
borrowed exclusively and retain their identities. Do not advance or mutate that
inner state outside this wrapper. Seed the fields then call `init!`; reinit
resets time, count and budget using the current fields as a new initial state.
It validates finite homogeneous PEC data without clamping or projecting them.

## Profile and boundary rules

`cells` gives nonnegative integer lower/upper thicknesses. Their sum must be less
than the axis cell count, retaining at least one interior cell. A positive width
requires positive finite peak decay rate after grid-T conversion; a disabled face
requires resolved rate zero. Periodic axes require both faces disabled. A fully
zero profile is allowed. Profiles are `:quadratic` and `:cubic`; no arbitrary
callback or inferred target reflection is supplied.

In each axis the profile rises from zero at the inner interface to its peak at
the outer wall. Node coordinate is i-1 and center coordinate is i-1/2, in cell
units. Profile fractions and their powers retain exact integer geometry. Each
E/B component selects node/center rates through the canonical Yee staggering.
At edges/corners the active axis rates add **before** forming one damping factor;
separate face masks are not multiplied and loss is counted once.

Outer boundaries remain the actual grid's PEC/periodic boundaries. Tangential E
and normal B remain zero on PEC faces/intersections, because damping preserves
zero and the canonical curl map preserves the same constraints. Periodic axes
have no duplicate endpoint and no layer in that direction. The method does not
turn a conducting wall into an exact open boundary.

## Equations, map and admission

With epsilon=1 and mu=1/c², the chosen artificial equations are

    dE/dt = c² curl(B) - lambda E - J_external
    dB/dt = -curl(E) - lambda B.

`peak_rate` specifies lambda in inverse time units: electric conductivity is
lambda and magnetic conductivity is lambda/c². The code applies equal decay
rates directly, without forming that quotient. All field components are damped.

Each fixed positive converted h uses conductivity midpoint half-steps around
the unchanged vacuum candidate. For each half-step, r=h*lambda/4 and
q=(1-r)/(1+r); output is the shared nearest/even rounding of q times the original
stored field. Original full h, profile and factor arithmetic remain rational
until that final output. A separately rounded half-duration, factor or product
cannot erase a representable completed result. Zero samples retain their signs.

Require h below the existing strict cached Maxwell CFL and the conservative
exact bound `h*sum(max(lower_peak,upper_peak))/4 < 1`. The latter excludes the
midpoint sign-flip regime and guarantees 0<q<=1 at all samples. Nearest rounding
and gradual underflow are required and never changed. This local damping bound
is not a global stability theorem for the split field map. Variable h is not
supported. With all rates zero, the identical canonical vacuum candidate acts
directly on the original fields; no damping arithmetic is performed.

`step_maxwell!(layer,h,impulse)` retains existing D1-only prescribed-current
support, shapes, precision and PEC mask. Its impulse arrays must be mutually
independent and independent of all 18 physical/candidate/predamped field arrays
and derived axis caches; this extends the vacuum model's storage admission.
Its work uses the fields
immediately before/after the middle sourced KDK map. The source can overlap the
layer, but there is no inferred source waveform accuracy, particle injection or
consistent charge-current model. D2/3 prescribed impulses reject.

## Energy and divergence accounting

Use the canonical physical metric H=prod(stored_dx)*sum(E²+c²B²)/2. Let F0 be
the old fields, Fa pre-damped, Fb after the canonical middle map, and F1 after
post-damping. Exact stored-state diagnostics retain

    L = H(F0)-H(Fa)+H(Fb)-H(F1)
    R = H(Fb)-H(Fa)+W
    H(F1)-H(F0)+L+W = R.

L is actual energy removed by damping, nonnegative under the specified rounded
submap. W is canonical prescribed-source work, positive out of the fields.
R is the signed physical Yee-map balance defect, including its distinction from
modified energy and staging roundoff. It is generally nonzero. A zero exact
bookkeeping residual does not imply exact physical energy conservation or a
contractive full step. Loss is not outward Poynting flux: PEC flux remains zero,
and paired periodic flux cancels in the existing physical flux diagnostic.

`absorber_diagnostics` returns canonical rounded field energy, its exact rational
value, initial energy, last/cumulative loss, source work and field defect, clock
and bookkeeping residual. Rational ledgers preserve complete signed differences
even when rounded physical energy is Inf. Existing `em_field_energy`,
`yee_poynting` and `yee_poynting_flux` delegate their actual vacuum meanings.

Diagonal conductivity changes electric and magnetic divergence in the layer;
compatible curls alone still annihilate divergence. These artificial sink effects
are not projected away. This layer cannot be attached to coupled PIC unchanged:
particle/surface charge, magnetic-divergence and injection laws need a separate
derivation. Finite-grid reflection across frequency/angle/resolution and layer
thickness remains unmeasured, including face/edge/corner configurations.

## Transactions and restart

The model owns one extra predamped YeeFields and derived O(sum(n)) exact axis
caches. No full rational field matrices or PML memory are invented. Fields,
clock, phase and all exact ledger changes stage before shared fixed-size physical
publication. Failure changes only private scratch. Profile/cache/storage/source
independence and current exact budget identity validate before staging. Caches
validate against the same canonical profile owner; init/restore builds fresh ones.
The inner CFL cache must equal the canonical grid/speed-derived value.
Expected-error publication assumes ordinary storage and exclusive access;
resource exhaustion, protected memory, custom setters and concurrent mutation
are outside that guarantee. No allocation or speed result has been measured.

Distinct `:yee_conductive_absorber` schema1 retains the exact scalar-midpoint/
staggered-profile policies, fixed configuration, nested synchronized vacuum
schema1 state, ready/count metadata and all exposed exact budget ledgers. Incoming
types, fields, configuration and ledger identity validate/detach before nested
restore. All derived caches prepare before publication; old destination field,
clock and ledger values may be stale. Fields restore verbatim with no projection.
Only initialized ready states save. PML, particle and old vacuum-only payloads
cannot substitute. The checksummed archive container remains version4, and no
forcing waveform or RNG is invented. No tests were added or run; runtime/restart
and broader CAP-002C/D acceptance remain open.
