# integrator.jl — explicit hybrid leapfrog with a subcycled magnetic field.
#
# Documented time levels (the cycle below advances all of these by one dt):
#   particle position x        : integer n          (→ n+1)
#   particle velocity v        : half    n−1/2      (→ n+1/2)
#   magnetic field   B         : integer n          (→ n+1)
#   evolving electron pe        : integer n          (→ n+1, optional closure)
#   electric field   E (carried): integer n         (recomputed → n+1)
#   transport ion force E_force : integer n         (recomputed with E, TransportElectrons)
#   moments for the B-subcycle : n+1/2 (midpoint positions x^{n+1/2}, v^{n+1/2})
#
# The loaded velocity is the physical v^0, so the FIRST step primes the leapfrog once
# (v^{-1/2} = v^0 − (dt/2)·a^0) — without it the loaded v^0 acts as v^{-1/2} and the run is only
# 1st-order in dt. This scheme is 2nd-order in dt (verified by a convergence study, incl. on
# magnetized problems); getting there needs both the priming and the step-4 u_i re-centering below.
#
# One step (CAM/CL-style; the B subcycle is RK4 for stability of the stiff
# whistler branch, with the n+1/2 ion moments held frozen across the subcycle):
#   1. Gather E^n, B^n to particles; Boris push  v^{n−1/2}→v^{n+1/2}, x^n→x^{n+1}.
#      Record midpoint positions x^{n+1/2}; wrap positions (periodic box).
#   2. Deposit frozen moments n^{n+1/2}, u_i^{n+1/2} from (x^{n+1/2}, v^{n+1/2}).
#   3. Subcycle B: integrate ∂B/∂t = −∇×E(B; frozen moments) over dt with N_B RK4
#      substeps  →  B^{n+1}.
#   4. Recompute carried E from n^{n+1}, u_i^{n+1/2}, B^{n+1}; then RE-CENTER u_i to n+1 via a
#      predictor half-kick and recompute → 2nd-order E^{n+1} (u_i^{n+1/2} alone would be 1st-order
#      in the −u_i×B convection term). See _recenter_carried_E!.
# The Boris update is centered at n; the B update is centered at n+1/2; the scheme is 2nd-order.

# Evolving pressure alone owns these grids and rollback buffers. Algebraic
# closures use `nothing`, so their hot path carries no extra array storage.
# `Eforce` is the transport-owned carried ion force and `saved_Eforce` its
# rollback copy; both are `nothing` (zero storage) for EvolvingElectrons.
mutable struct _PressureEvolutionWorkspace{D,T,F}
    ue::NTuple{3,Array{T,D}}
    divue::Array{T,D}
    kp::NTuple{4,Array{T,D}}
    ptmp::Array{T,D}
    saved::HybridFields{D,T}
    x0::NTuple{D,Vector{T}}
    v0::NTuple{3,Vector{T}}
    clock::Tuple{T,Int,T}
    Eforce::F
    saved_Eforce::F
end

_transport_force_storage(::FourierGrid, ::ElectronClosure) = nothing
_transport_force_storage(g::FourierGrid{D,T}, ::TransportElectrons) where {D,T} =
    ntuple(_ -> zeros(T, g.n), 3)

function _pressure_workspace(g::FourierGrid{D,T}, closure, count) where {D,T}
    _evolves_pressure(closure) || return nothing
    force = _transport_force_storage(g, closure)
    return _PressureEvolutionWorkspace{D,T,typeof(force)}(
        ntuple(_ -> zeros(T, g.n), 3),
        zeros(T, g.n),
        ntuple(_ -> zeros(T, g.n), 4),
        zeros(T, g.n),
        HybridFields{D,T}(g.n),
        ntuple(_ -> zeros(T, count), D),
        ntuple(_ -> zeros(T, count), 3),
        (zero(T), 0, zero(T)),
        force,
        _transport_force_storage(g, closure),
    )
end


"""
    HybridStepper(g, model, shape, N::Integer)

Workspace and configuration for time-stepping ions on grid `g` under `model`
with deposition/gather `shape`. `N` is the initial total particle workspace
length, resized as needed. Pass one proton `ParticleSet` or a vector of positive
ion populations to [`init!`](@ref) and [`step!`](@ref), after setting the initial
field in `stepper.fields.B`. With [`EvolvingElectrons`](@ref), also set the
finite, nonnegative initial pressure in `stepper.fields.pe`; initialization
preserves that pressure. Other closures compute pressure from their algebraic law.
[`TransportElectrons`](@ref) uses the same pressure state and additionally
allocates the carried ion force and its rollback copy (six grid arrays).
"""
mutable struct HybridStepper{D,T,SH<:ShapeFunction,M<:HybridModel,W}
    g::FourierGrid{D,T}
    model::M
    shape::SH
    fields::HybridFields{D,T}          # canonical B, carried E, plus Ohm scratch
    fn::Array{T,D}                     # frozen n^{n+1/2}
    fui::NTuple{3,Array{T,D}}          # frozen u_i^{n+1/2}
    Escr::NTuple{3,Array{T,D}}         # E scratch for the B-subcycle
    k1::NTuple{3,Array{T,D}}
    k2::NTuple{3,Array{T,D}}
    k3::NTuple{3,Array{T,D}}
    k4::NTuple{3,Array{T,D}}
    Btmp::NTuple{3,Array{T,D}}
    Ep::NTuple{3,Vector{T}}
    Bp::NTuple{3,Vector{T}}
    vpred::NTuple{3,Vector{T}}         # predicted v^{n+1} for the carried-E re-centering (2nd order)
    xmid::NTuple{D,Vector{T}}
    work::Vector{T}
    time::Base.RefValue{T}
    step::Base.RefValue{Int}
    last_dt::Base.RefValue{T}          # previous drift interval fixes the velocity half level
    pressure_work::W
end

const _HybridPopulations{D,T} = Union{ParticleSet{D,T},AbstractVector{<:ParticleSet{D,T}}}

function HybridStepper(
    g::FourierGrid{D,T},
    model::M,
    shape::SH,
    N::Integer,
) where {D,T,M<:HybridModel,SH<:ShapeFunction}
    nc = g.n
    Np = _particle_length(N)
    _validate_evolving_model(model, T)
    pressure_work = _pressure_workspace(g, model.closure, Np)
    HybridStepper{D,T,SH,M,typeof(pressure_work)}(
        g,
        model,
        shape,
        HybridFields{D,T}(nc; anisotropic = is_anisotropic(model.closure)),
        zeros(T, nc),
        ntuple(_ -> zeros(T, nc), Val(3)),
        ntuple(_ -> zeros(T, nc), Val(3)),
        ntuple(_ -> zeros(T, nc), Val(3)),
        ntuple(_ -> zeros(T, nc), Val(3)),
        ntuple(_ -> zeros(T, nc), Val(3)),
        ntuple(_ -> zeros(T, nc), Val(3)),
        ntuple(_ -> zeros(T, nc), Val(3)),
        ntuple(_ -> zeros(T, Np), Val(3)),   # Ep
        ntuple(_ -> zeros(T, Np), Val(3)),   # Bp
        ntuple(_ -> zeros(T, Np), Val(3)),   # vpred
        ntuple(_ -> zeros(T, Np), Val(D)),   # xmid
        zeros(T, Np),
        Ref(zero(T)),
        Ref(0),
        Ref(zero(T)),                        # last_dt
        pressure_work,
    )
end

_pressure_populations(ps::ParticleSet) = (ps,)
_pressure_populations(species::AbstractVector) = species

function _validate_pressure_storage(st::HybridStepper, ps)
    w = st.pressure_work
    w === nothing && return nothing
    f, g = st.fields, st.g
    _validate_evolving_pressure_storage(f, g)
    for a in (
        st.fn,
        st.fui...,
        st.Escr...,
        st.k1...,
        st.k2...,
        st.k3...,
        st.k4...,
        st.Btmp...,
        w.ue...,
        w.divue,
        w.kp...,
        w.ptmp,
        w.saved.n,
        w.saved.ui...,
        w.saved.B...,
        w.saved.E...,
        w.saved.J...,
        w.saved.pe,
        w.saved.gradp...,
        w.saved.ninv,
        w.saved.lapJ...,
        st.Ep...,
        st.Bp...,
        st.vpred...,
        st.xmid...,
        st.work,
        w.x0...,
        w.v0...,
    )
        _particle_arrays_mightalias(f.pe, a) &&
            throw(ArgumentError("evolving pressure must have independent storage"))
    end
    for population in _pressure_populations(ps),
        a in (population.x..., population.v..., population.weight, population.id, population.tag)

        _particle_arrays_mightalias(f.pe, a) &&
            throw(ArgumentError("evolving pressure must not alias particle storage"))
    end
    _validate_transport_force_storage(st, ps, w.Eforce)
    return nothing
end

_validate_transport_force_storage(st, ps, ::Nothing) = nothing

# The carried force and its rollback copy must be independent grid arrays, as
# for stored pressure: aliasing would corrupt the paired field or rollback.
# Grid and particle-length storage are scanned as separate homogeneous tuples,
# so the per-step check stays type-stable and allocation-free in every D.
function _validate_transport_force_storage(
    st::HybridStepper{D,T},
    ps,
    force::NTuple{3,Array{T,D}},
) where {D,T}
    w, f = st.pressure_work, st.fields
    carried = (force..., w.saved_Eforce...)
    # The field type fixes Array{T,D} with one-based axes, so matching the
    # already grid-validated pressure size is the complete layout check.
    for a in carried
        size(a) == size(f.pe) || throw(
            DimensionMismatch(
                "transport force size $(size(a)) does not match grid size $(size(f.pe))",
            ),
        )
    end
    for (i, a) in enumerate(carried), j = 1:(i-1)
        _particle_arrays_mightalias(a, carried[j]) &&
            throw(ArgumentError("transport force arrays must have independent storage"))
    end
    grids = (
        f.n,
        f.ui...,
        f.B...,
        f.E...,
        f.J...,
        f.pe,
        f.gradp...,
        f.ninv,
        f.lapJ...,
        f.pforce...,
        f.cgl_dp,
        st.fn,
        st.fui...,
        st.Escr...,
        st.k1...,
        st.k2...,
        st.k3...,
        st.k4...,
        st.Btmp...,
        w.ue...,
        w.divue,
        w.kp...,
        w.ptmp,
        w.saved.n,
        w.saved.ui...,
        w.saved.B...,
        w.saved.E...,
        w.saved.J...,
        w.saved.pe,
        w.saved.gradp...,
        w.saved.ninv,
        w.saved.lapJ...,
    )
    scratch = (st.Ep..., st.Bp..., st.vpred..., st.xmid..., st.work, w.x0..., w.v0...)
    for b in grids, a in carried
        _particle_arrays_mightalias(a, b) &&
            throw(ArgumentError("transport force must have independent storage"))
    end
    for b in scratch, a in carried
        _particle_arrays_mightalias(a, b) &&
            throw(ArgumentError("transport force must have independent storage"))
    end
    for population in _pressure_populations(ps), a in carried
        aliased =
            any(b -> _particle_arrays_mightalias(a, b), population.x) ||
            any(b -> _particle_arrays_mightalias(a, b), population.v) ||
            _particle_arrays_mightalias(a, population.weight) ||
            _particle_arrays_mightalias(a, population.id) ||
            _particle_arrays_mightalias(a, population.tag)
        aliased && throw(ArgumentError("transport force must not alias particle storage"))
    end
    return nothing
end

function _validate_pressure_step_state(st::HybridStepper{D,T}, ps) where {D,T}
    st.pressure_work === nothing && return nothing
    _validate_evolving_model(st.model, T)
    _validate_pressure_storage(st, ps)
    # The transaction copies x/v before deposition; validate their canonical
    # layout here so MPI ranks reject malformed storage collectively.
    for population in _pressure_populations(ps)
        _validate_moment_particles(population)
    end
    _validate_evolving_pressure(st.fields.pe, st.g)
    for a in (st.fields.B..., st.fields.E...)
        _validate_pressure_finite(a, "evolving-pressure field input")
    end
    _validate_transport_step_inputs(st, ps, st.pressure_work.Eforce)
    return nothing
end

_validate_transport_step_inputs(st, ps, ::Nothing) = nothing

function _validate_transport_step_inputs(st, ps, force::NTuple{3})
    # The scalar interface deposits proton-normalized number density; only a
    # unit charge makes it the charge density the paired energy balance needs.
    ps isa ParticleSet &&
        !isone(ps.q) &&
        throw(
            ArgumentError(
                "TransportElectrons scalar populations require q = 1; use a one-element collection",
            ),
        )
    for a in force
        _validate_pressure_finite(a, "transport carried ion force")
    end
    return nothing
end

_copy_transport_force!(::Nothing, ::Nothing) = nothing
function _copy_transport_force!(destination::NTuple{3}, source::NTuple{3})
    foreach(copyto!, destination, source)
    return nothing
end

function _save_pressure_transaction!(st::HybridStepper{D,T}, ps) where {D,T}
    w = st.pressure_work
    w === nothing && return nothing
    total_count = sum(nparticles(p) for p in _pressure_populations(ps); init = 0)
    if length(w.v0[1]) != total_count
        w.x0 = ntuple(_ -> Vector{T}(undef, total_count), D)
        w.v0 = ntuple(_ -> Vector{T}(undef, total_count), 3)
    end
    _copy_hybrid_fields!(w.saved, st.fields)
    _copy_transport_force!(w.saved_Eforce, w.Eforce)
    offset = 0
    for population in _pressure_populations(ps)
        population_count = nparticles(population)
        for d = 1:D
            copyto!(w.x0[d], offset + 1, population.x[d], 1, population_count)
        end
        for c = 1:3
            copyto!(w.v0[c], offset + 1, population.v[c], 1, population_count)
        end
        offset += population_count
    end
    w.clock = (st.time[], st.step[], st.last_dt[])
    return nothing
end

function _restore_pressure_transaction!(st::HybridStepper{D,T}, ps) where {D,T}
    w = st.pressure_work
    w === nothing && return nothing
    _copy_hybrid_fields!(st.fields, w.saved)
    _copy_transport_force!(w.Eforce, w.saved_Eforce)
    offset = 0
    for population in _pressure_populations(ps)
        count = nparticles(population)
        for d = 1:D
            copyto!(population.x[d], 1, w.x0[d], offset + 1, count)
        end
        for c = 1:3
            copyto!(population.v[c], 1, w.v0[c], offset + 1, count)
        end
        offset += count
    end
    st.time[], st.step[], st.last_dt[] = w.clock
    return nothing
end

function _validate_pressure_electric!(st::HybridStepper)
    st.pressure_work === nothing && return nothing
    for a in st.fields.E
        _validate_pressure_finite(a, "evolving-pressure carried electric field")
    end
    _validate_transport_force_finite(st.pressure_work.Eforce)
    return nothing
end

_validate_transport_force_finite(::Nothing) = nothing
function _validate_transport_force_finite(force::NTuple{3})
    for a in force
        _validate_pressure_finite(a, "transport carried ion force")
    end
    return nothing
end

function _resize_hybrid_particle_workspaces!(st, n::Integer)
    n >= 0 || throw(ArgumentError("particle workspace length must be nonnegative, got $n"))
    N = Int(n)
    D = length(st.xmid)
    T = eltype(st.work)
    length(st.work) == N &&
        all(length(st.Ep[c]) == N for c = 1:3) &&
        all(length(st.Bp[c]) == N for c = 1:3) &&
        all(length(st.vpred[c]) == N for c = 1:3) &&
        all(length(st.xmid[d]) == N for d = 1:D) &&
        return st

    st.Ep = ntuple(_ -> Vector{T}(undef, N), 3)
    st.Bp = ntuple(_ -> Vector{T}(undef, N), 3)
    st.vpred = ntuple(_ -> Vector{T}(undef, N), 3)
    st.xmid = ntuple(_ -> Vector{T}(undef, N), D)
    st.work = Vector{T}(undef, N)
    return st
end

function _validate_hybrid_species(ps::ParticleSet)
    _validated_qm(ps)
    return nparticles(ps)
end

function _require_independent_hybrid_species(species)
    for (i, ps) in enumerate(species)
        _particle_components_mightalias(ps) &&
            throw(ArgumentError("hybrid particle components must not alias"))
        for (j, other) in enumerate(species)
            j >= i && break
            _particle_sets_mightalias(ps, other) &&
                throw(ArgumentError("hybrid populations must have independent particle storage"))
        end
    end
    return nothing
end

function _validate_hybrid_species(species::AbstractVector{<:ParticleSet{D,T}}) where {D,T}
    _require_independent_hybrid_species(species)
    total = 0
    for ps in species
        _validated_qm(ps)
        ps.q > zero(T) || throw(ArgumentError("hybrid ion charge must be positive"))
        n = nparticles(ps)
        for a in (ps.x..., ps.v..., ps.weight, ps.id, ps.tag)
            length(a) == n || throw(DimensionMismatch("hybrid particle arrays must have length $n"))
            _check_particle_vector_axes(:hybrid_particle, a, n)
        end
        _validate_particle_positions(ps)
        _require_finite_particle_velocities(ps, "hybrid stepping")
        for w in ps.weight
            isfinite(w) && w >= zero(T) ||
                throw(ArgumentError("hybrid particle weights must be finite and nonnegative"))
        end
        n <= typemax(Int) - total || throw(ArgumentError("total particle count exceeds Int"))
        total += n
    end
    return total
end

# Each population gets a contiguous part of the existing total-length buffers.
function _hybrid_particle_workspace(st, indices)
    return (
        Ep = map(a -> view(a, indices), st.Ep),
        Bp = map(a -> view(a, indices), st.Bp),
        vpred = map(a -> view(a, indices), st.vpred),
        xmid = map(a -> view(a, indices), st.xmid),
        work = view(st.work, indices),
    )
end

# Serial and MPI callers share raw moments; MPI reduces across ranks before
# either caller converts the momentum density to bulk velocity.
function _hybrid_integer_raw_moments!(st, ps::ParticleSet{D,T}) where {D,T}
    _deposit_density_momentum_moments!(st.fields.n, st.fields.ui, ps, st.g, st.shape, st.work)
    return st
end

function _hybrid_integer_raw_moments!(st, species::AbstractVector{<:ParticleSet{D,T}}) where {D,T}
    f = st.fields
    _zero_ion_moments!(f.n, f.ui)
    offset = 0
    for ps in species
        n = nparticles(ps)
        work = view(st.work, (offset+1):(offset+n))
        _accumulate_charge_moments!(f.n, f.ui, ps, st.g, st.shape, st.fn, st.fui, work)
        offset += n
    end
    return st
end

function _hybrid_integer_moments!(st, ps)
    _hybrid_integer_raw_moments!(st, ps)
    _normalize_ion_moments!(st.fields.n, st.fields.ui, eltype(st.fields.n)(st.model.nfloor))
    return st
end

"""
    init!(stepper, ps)
    init!(stepper, species::AbstractVector)

Compute the carried electric field E⁰ from the initial particle moments and
`stepper.fields.B`. Call once after loading particles and setting B. With
[`EvolvingElectrons`](@ref), set finite nonnegative `stepper.fields.pe` first;
initialization preserves it, including an exact zero cold state. Reinitialization
also preserves pressure and resets the clock and leapfrog phase. A rejected
evolving-pressure initialization restores public fields and clocks.
Collection calls combine positive ion populations with their own charge and mass
using charge-weighted moments. Populations require independent storage and a
common time level. Reinitialize when adding populations loaded at a new physical
velocity time level. Empty populations contribute zero moments.
[`TransportElectrons`](@ref) also computes the carried ion force and requires
the initial density to be finite and at least `nfloor` in every cell.
"""
function init!(st::HybridStepper{D,T}, ps::_HybridPopulations{D,T}) where {D,T}
    count = ps isa ParticleSet ? nparticles(ps) : _validate_hybrid_species(ps)
    _validate_pressure_step_state(st, ps)
    _save_pressure_transaction!(st, ps)
    completed = false
    try
        _resize_hybrid_particle_workspaces!(st, count)
        _hybrid_integer_moments!(st, ps)
        _carried_ohms_law!(st)
        _validate_pressure_electric!(st)
        st.time[] = zero(T)                 # (re)start at t=0; step==0 triggers leapfrog priming
        st.step[] = 0
        st.last_dt[] = zero(T)
        completed = true
        return st
    finally
        # Cleanup only: the original exception/stack propagates, including interrupts.
        completed || _restore_pressure_transaction!(st, ps)
    end
end

# F = −∇×E(Btrial; frozen moments) into out
@inline function _bfield_rhs!(
    out::NTuple{3,<:Array{T,D}},
    Btrial::NTuple{3,<:Array{T,D}},
    st::HybridStepper{D,T},
) where {D,T}
    f = st.fields
    # frozen pe/∇pe/ninv (scalar) or the frozen ∇·P_e force (anisotropic/CGL) were
    # precomputed once per step; only the B-dependent part is evaluated here each subcycle
    # stage. The `is_anisotropic` test is compile-time-resolved from the closure type, so
    # scalar closures keep the exact original call below.
    if is_anisotropic(st.model.closure)
        # the gyrotropic force depends on the field DIRECTION b, so it must track the trial
        # B within the subcycle (a step-frozen b would misalign as B_⊥ grows and drive a
        # numerical instability); recompute ∇·P_e(n^{n+1/2}, B_trial) each stage.
        _anisotropic_pressure_force!(f, st.fn, Btrial, st.model.closure, st.g)
        _ohm_Efield_aniso!(
            st.Escr,
            st.fui,
            Btrial,
            f.J,
            f.lapJ,
            f.pforce,
            f.ninv,
            T(st.model.η),
            T(st.model.ηH),
            st.g,
        )
    else
        _ohm_Efield!(
            st.Escr,
            st.fui,
            Btrial,
            f.J,
            f.lapJ,
            f.gradp,
            f.ninv,
            T(st.model.η),
            T(st.model.ηH),
            st.g,
        )
    end
    _apply_electron_inertia!(st.Escr, T(st.model.de2), st.g)   # E ← E/(1+d_e²k²) if de2>0
    faraday_rhs!(out, st.Escr, st.g)
    return out
end

# one RK4 substep of size h on canonical B (st.fields.B), frozen moments
function _rk4_B!(st::HybridStepper{D,T}, h::T) where {D,T}
    B = st.fields.B
    _bfield_rhs!(st.k1, B, st)
    for c = 1:3
        @. st.Btmp[c] = B[c] + (h / 2) * st.k1[c]
    end
    _bfield_rhs!(st.k2, st.Btmp, st)
    for c = 1:3
        @. st.Btmp[c] = B[c] + (h / 2) * st.k2[c]
    end
    _bfield_rhs!(st.k3, st.Btmp, st)
    for c = 1:3
        @. st.Btmp[c] = B[c] + h * st.k3[c]
    end
    _bfield_rhs!(st.k4, st.Btmp, st)
    for c = 1:3
        @. B[c] += (h / 6) * (st.k1[c] + 2 * st.k2[c] + 2 * st.k3[c] + st.k4[c])
    end
    return B
end

# Constant Ohmic transport adds ηJ to the stage field that drives Faraday and
# its heat (γ−1)η|J|² to the pressure rate, both from this stage's current.
# The ideal closure and η = 0 leave the stage arithmetic untouched.
_add_transport_stage_sources!(E, rate, J, model::HybridModel) =
    _add_transport_stage_sources!(E, rate, J, model, model.closure)
_add_transport_stage_sources!(E, rate, J, model, ::ElectronClosure) = nothing
function _add_transport_stage_sources!(
    E,
    rate::Array{T},
    J,
    model,
    closure::TransportElectrons,
) where {T}
    η = T(model.η)
    iszero(η) && return nothing
    _add_ohmic_field!(E, E, J, η)
    _add_ohmic_pressure_heating!(rate, J, T(closure.γ), η)
    return nothing
end

# One coupled field RHS: E/J and pressure transport use the same B/p stage.
function _bpressure_rhs!(
    outB,
    outp,
    Btrial,
    ptrial,
    st::HybridStepper{D,T},
    cycle,
    stage,
) where {D,T}
    f, w, g = st.fields, st.pressure_work, st.g
    _validate_evolving_pressure(ptrial, g; subcycle = cycle, stage)
    for a in Btrial
        _validate_pressure_finite(a, "evolving-pressure magnetic stage")
    end
    gradient!(f.gradp, ptrial, g)
    _ohm_Efield!(st.Escr, st.fui, Btrial, f.J, f.lapJ, f.gradp, f.ninv, zero(T), zero(T), g)
    electron_velocity!(w.ue, st.fui, f.J, st.fn; nfloor = T(st.model.nfloor))
    divergence!(w.divue, w.ue, g)
    _electron_pressure_rhs!(outp, ptrial, w.ue, f.gradp, w.divue, T(st.model.closure.γ), Val(D))
    _add_transport_stage_sources!(st.Escr, outp, f.J, st.model)
    faraday_rhs!(outB, st.Escr, g)
    _validate_pressure_finite(outp, "evolving-pressure rate")
    for a in outB
        _validate_pressure_finite(a, "evolving-pressure magnetic rate")
    end
    return nothing
end

function _rk4_Bp!(st::HybridStepper{D,T}, h::T, cycle) where {D,T}
    B, pe, w = st.fields.B, st.fields.pe, st.pressure_work
    _bpressure_rhs!(st.k1, w.kp[1], B, pe, st, cycle, 1)
    for c = 1:3
        @. st.Btmp[c] = B[c] + (h / 2) * st.k1[c]
    end
    @. w.ptmp = pe + (h / 2) * w.kp[1]
    _bpressure_rhs!(st.k2, w.kp[2], st.Btmp, w.ptmp, st, cycle, 2)
    for c = 1:3
        @. st.Btmp[c] = B[c] + (h / 2) * st.k2[c]
    end
    @. w.ptmp = pe + (h / 2) * w.kp[2]
    _bpressure_rhs!(st.k3, w.kp[3], st.Btmp, w.ptmp, st, cycle, 3)
    for c = 1:3
        @. st.Btmp[c] = B[c] + h * st.k3[c]
    end
    @. w.ptmp = pe + h * w.kp[3]
    _bpressure_rhs!(st.k4, w.kp[4], st.Btmp, w.ptmp, st, cycle, 4)
    for c = 1:3
        @. st.Btmp[c] = B[c] + (h / 6) * (st.k1[c] + 2 * st.k2[c] + 2 * st.k3[c] + st.k4[c])
        _validate_pressure_finite(st.Btmp[c], "evolving-pressure magnetic candidate")
    end
    @. w.ptmp = pe + (h / 6) * (w.kp[1] + 2 * w.kp[2] + 2 * w.kp[3] + w.kp[4])
    _validate_evolving_pressure(w.ptmp, st.g; subcycle = cycle, stage = 5)
    for c = 1:3
        copyto!(B[c], st.Btmp[c])
    end
    copyto!(pe, w.ptmp)
    return nothing
end

# Serial and replicated MPI share preparation and every subcycle. Algebraic
# closures retain their original pressure preparation and B-only RK arithmetic.
function _advance_hybrid_fields!(st::HybridStepper{D,T}, dtT::T, NB::Integer) where {D,T}
    f, g = st.fields, st.g
    if _evolves_pressure(st.model.closure)
        _validate_pressure_finite(st.fn, "evolving-pressure midpoint density")
        for u in st.fui
            _validate_pressure_finite(u, "evolving-pressure midpoint ion velocity")
        end
        _transports_electrons(st.model.closure) &&
            _validate_transport_density(st.fn, T(st.model.nfloor), "midpoint")
        _ohm_ninv!(f.ninv, st.fn, T(st.model.nfloor), f.floor_count)
    elseif is_anisotropic(st.model.closure)
        _ohm_ninv!(f.ninv, st.fn, T(st.model.nfloor), f.floor_count)
    else
        _ohm_prep!(
            f.pe,
            f.gradp,
            f.ninv,
            st.fn,
            st.model.closure,
            T(st.model.nfloor),
            f.floor_count,
            g,
        )
    end
    hb = dtT / NB
    for cycle = 1:NB
        if _evolves_pressure(st.model.closure)
            _rk4_Bp!(st, hb, cycle)
        else
            _rk4_B!(st, hb)
        end
    end
    return nothing
end

function _validated_step_dt(
    ::Type{T},
    dt::Real,
    NB::Integer;
    min_NB::Integer,
    name::AbstractString,
) where {T}
    NB >= min_NB ||
        throw(ArgumentError("$name requires NB >= $min_NB magnetic subcycles (got $NB)"))
    return _validated_nonnegative_dt(T, dt; name)
end

# One-time leapfrog priming. The loaded velocity is the physical v^0, but the leapfrog Boris push
# expects v^{-1/2}; push v back a half step v^{-1/2} = v^0 − (dt/2)·a^0 (forward-Euler, with
# a^0 = qm(E^0 + v^0×B^0) from the just-gathered fields) so the first Boris kick is centred at
# t=0. Applied once (first step after init!) — a one-time O(dt²) IC correction that restores the
# documented 2nd-order accuracy for a real init!+step! run (verified: end-to-end rate 1 → 2).
# `rel=true` (relativistic push, works in momentum u=γv) backs up MOMENTUM: u^{-1/2}=γ^0 v^0 − h·a^0,
# v^{-1/2}=u^{-1/2}/γ(u); the default non-relativistic path is byte-identical for the hybrid callers.
@inline function _prime_leapfrog!(v, Ep, Bp, qm, h, np, rel::Bool = false, c = one(eltype(v[1])))
    vx, vy, vz = v
    ex, ey, ez = Ep
    bx, by, bz = Bp
    if rel
        c2 = c * c
        o = one(c)
        @inbounds for p = 1:np
            v0x, v0y, v0z = vx[p], vy[p], vz[p]
            βx = v0x / c
            βy = v0y / c
            βz = v0z / c
            β2 = βx * βx + βy * βy + βz * βz
            (isfinite(β2) && β2 < o) || throw(
                ArgumentError(
                    "_prime_leapfrog!: input velocity must be finite and satisfy |v| < c",
                ),
            )
            g0 = o / sqrt(o - β2)
            ax = qm * (ex[p] + (v0y * bz[p] - v0z * by[p]))
            ay = qm * (ey[p] + (v0z * bx[p] - v0x * bz[p]))
            az = qm * (ez[p] + (v0x * by[p] - v0y * bx[p]))
            umx = g0 * v0x - h * ax
            umy = g0 * v0y - h * ay
            umz = g0 * v0z - h * az
            gm = sqrt(o + (umx * umx + umy * umy + umz * umz) / c2)
            vx[p] = umx / gm
            vy[p] = umy / gm
            vz[p] = umz / gm
        end
    else
        @inbounds for p = 1:np
            ux, uy, uz = vx[p], vy[p], vz[p]
            vx[p] = ux - h * qm * (ex[p] + (uy * bz[p] - uz * by[p]))
            vy[p] = uy - h * qm * (ey[p] + (uz * bx[p] - ux * bz[p]))
            vz[p] = uz - h * qm * (ez[p] + (ux * by[p] - uy * bx[p]))
        end
    end
    return v
end

# Carried fields. The full E drives Faraday's law and diagnostics. With
# TransportElectrons the ions are pushed by the carried force E_force instead:
# the first-step priming, every Boris kick and the endpoint predictor gather it.
# Only HybridStepper owns that force; CAM-CL rejects the closure on construction.
_particle_electric_field(st) = _particle_electric_field(st, st.model.closure)
_particle_electric_field(st, ::ElectronClosure) = st.fields.E
_particle_electric_field(st::HybridStepper, ::TransportElectrons) = st.pressure_work.Eforce

_carried_ohms_law!(st) = _carried_ohms_law!(st, st.model.closure)
_carried_ohms_law!(st, ::ElectronClosure) = ohms_law!(st.fields, st.model, st.g)
_carried_ohms_law!(st::HybridStepper, ::TransportElectrons) =
    _carried_transport_ohms_law!(st.fields, st.pressure_work, st.model, st.g)
# The grid field is abstractly typed; passing the mutable workspace instead of
# its force tuple keeps that dynamic call boundary free of a boxed argument.
_carried_transport_ohms_law!(f, w, model, g) = _transport_ohms_law!(f, w.Eforce, model, g)

# Forward half-kick predictor v^{n+1/2} → v^{n+1} into `vp`, from particle-gathered E, B.
@inline function _predict_half_kick!(vp, v, Ep, Bp, qm, h, np)
    vx, vy, vz = v
    vpx, vpy, vpz = vp
    ex, ey, ez = Ep
    bx, by, bz = Bp
    @inbounds for p = 1:np
        ux, uy, uz = vx[p], vy[p], vz[p]
        vpx[p] = ux + h * qm * (ex[p] + (uy * bz[p] - uz * by[p]))
        vpy[p] = uy + h * qm * (ey[p] + (uz * bx[p] - ux * bz[p]))
        vpz[p] = uz + h * qm * (ez[p] + (ux * by[p] - uy * bx[p]))
    end
    return vp
end

# Re-center the carried E to integer level n+1. The E just computed in step 4 deposits u_i from
# the half-step velocity v^{n+1/2}, so the −u_i×B convection term lags by dt/2 — an O(dt) force
# error that makes the whole scheme 1st-order on magnetized problems. Predict v^{n+1} with a
# forward half-kick from the fresh E^{n+1}, B^{n+1}, re-deposit u_i^{n+1}, and recompute the
# carried E, restoring the documented CAM/CL 2nd-order accuracy (verified: magnetized temporal
# convergence rate 1 → 2). Unmagnetized behavior is unchanged (u_i drops out of Ohm's law when
# B=0). Generic over the stepper (HybridStepper/CAMCLStepper share Ep/Bp/vpred/work). No alloc.
function _hybrid_predictor_raw_moments!(st, ps::ParticleSet{D,T}, dtT::T) where {D,T}
    qm = _validated_qm(ps)
    g = st.g
    gather_vector!(st.Ep, _particle_electric_field(st), ps, g, st.shape)   # E^{n+1}
    gather_vector!(st.Bp, st.fields.B, ps, g, st.shape)   # B^{n+1}
    _predict_half_kick!(st.vpred, ps.v, st.Ep, st.Bp, qm, dtT / 2, length(ps.weight))
    pspred = ParticleSet{D,T}(ps.x, st.vpred, ps.weight, ps.id, ps.tag, ps.q, ps.m)
    _deposit_density_momentum_moments!(st.fields.n, st.fields.ui, pspred, g, st.shape, st.work)
    return st
end

function _hybrid_predictor_raw_moments!(
    st,
    species::AbstractVector{<:ParticleSet{D,T}},
    dtT::T,
) where {D,T}
    f = st.fields
    _zero_ion_moments!(f.n, f.ui)
    offset = 0
    for ps in species
        n = nparticles(ps)
        ws = _hybrid_particle_workspace(st, (offset+1):(offset+n))
        gather_vector!(ws.Ep, _particle_electric_field(st), ps, st.g, st.shape)
        gather_vector!(ws.Bp, f.B, ps, st.g, st.shape)
        _predict_half_kick!(ws.vpred, ps.v, ws.Ep, ws.Bp, _validated_qm(ps), dtT / 2, n)
        predicted = ParticleSet{D,T}(ps.x, ws.vpred, ps.weight, ps.id, ps.tag, ps.q, ps.m)
        _accumulate_charge_moments!(f.n, f.ui, predicted, st.g, st.shape, st.fn, st.fui, ws.work)
        offset += n
    end
    return st
end

function _recenter_carried_E!(st, ps, dtT)
    _hybrid_predictor_raw_moments!(st, ps, dtT)
    _normalize_ion_moments!(st.fields.n, st.fields.ui, eltype(st.fields.n)(st.model.nfloor))
    _carried_ohms_law!(st)
    return st
end

# Shared gather/prime/kick/drift kernel. Fields stay fixed across populations.
function _push_hybrid_population!(st, ws, ps::ParticleSet{D,T}, dtT::T, kick_dt::T) where {D,T}
    qm = _validated_qm(ps)
    h = dtT / 2
    gather_vector!(ws.Ep, _particle_electric_field(st), ps, st.g, st.shape)
    gather_vector!(ws.Bp, st.fields.B, ps, st.g, st.shape)
    st.step[] == 0 && _prime_leapfrog!(ps.v, ws.Ep, ws.Bp, qm, h, length(ps.weight))
    vx, vy, vz = ps.v
    @inbounds for p in eachindex(ps.weight)
        nx, ny, nz = boris_kick(
            vx[p],
            vy[p],
            vz[p],
            ws.Ep[1][p],
            ws.Ep[2][p],
            ws.Ep[3][p],
            ws.Bp[1][p],
            ws.Bp[2][p],
            ws.Bp[3][p],
            qm,
            kick_dt,
        )
        vx[p] = nx
        vy[p] = ny
        vz[p] = nz
        for d = 1:D
            ws.xmid[d][p] = ps.x[d][p] + h * ps.v[d][p]
            ps.x[d][p] += dtT * ps.v[d][p]
        end
    end
    lo = ntuple(_ -> zero(T), D)
    apply_periodic!(ps, lo, st.g.L)
    midpoint = ParticleSet{D,T}(ws.xmid, ps.v, ps.weight, ps.id, ps.tag, ps.q, ps.m)
    apply_periodic!(midpoint, lo, st.g.L)
    return midpoint
end

function _push_hybrid_midpoint_raw_moments!(
    st,
    ps::ParticleSet{D,T},
    dtT::T,
    kick_dt::T,
) where {D,T}
    midpoint = _push_hybrid_population!(st, st, ps, dtT, kick_dt)
    _deposit_density_momentum_moments!(st.fn, st.fui, midpoint, st.g, st.shape, st.work)
    return st
end

function _push_hybrid_midpoint_raw_moments!(
    st,
    species::AbstractVector{<:ParticleSet{D,T}},
    dtT::T,
    kick_dt::T,
) where {D,T}
    _zero_ion_moments!(st.fn, st.fui)
    offset = 0
    for ps in species
        n = nparticles(ps)
        ws = _hybrid_particle_workspace(st, (offset+1):(offset+n))
        midpoint = _push_hybrid_population!(st, ws, ps, dtT, kick_dt)
        _accumulate_charge_moments!(
            st.fn,
            st.fui,
            midpoint,
            st.g,
            st.shape,
            st.fields.n,
            st.fields.ui,
            ws.work,
        )
        offset += n
    end
    return st
end

function _push_hybrid_midpoint_moments!(st, ps, dtT, kick_dt)
    _push_hybrid_midpoint_raw_moments!(st, ps, dtT, kick_dt)
    _normalize_ion_moments!(st.fn, st.fui, eltype(st.fn)(st.model.nfloor))
    return st
end

"""
    step!(stepper, ps, dt; NB=1)
    step!(stepper, species::AbstractVector, dt; NB=1)

Advance the plasma one timestep `dt` with `NB` magnetic subcycles (periodic box).
The particle kick spans the previous and next half-time levels when `dt` changes;
positions, the magnetic field and evolving pressure advance by the requested
`dt` to integer time; particle velocity remains at its half-time level.
Collection calls push every ion population with its own charge-to-mass ratio;
one field update uses their combined charge-weighted midpoint moments.

With [`EvolvingElectrons`](@ref), pressure and B advance through the same RK4
stages with frozen midpoint ion moments. This closure requires an ideal massless
model. Pressure must remain finite and nonnegative at every stage; a rejected
step restores particles, public fields and clocks, preserving their array
identities. There is no pressure clipping. Choose a smaller timestep or improve
spatial resolution when a pressure stage is rejected. A zero timestep changes
nothing, including invalid input pressure. Private work arrays are scratch.

[`TransportElectrons`](@ref) shares these stages and adds `ηJ` to the stage
field and `(γ−1)η|J|²` to the stage pressure rate; ions are pushed by the
carried force `E_force`. The midpoint, endpoint and predictor densities must
stay finite and at least `nfloor`, and the transaction also restores the
carried force.
"""
function step!(
    st::HybridStepper{D,T},
    ps::_HybridPopulations{D,T},
    dt::Real;
    NB::Integer = 1,
) where {D,T}
    dtT = _validated_step_dt(T, dt, NB; min_NB = 1, name = "step!")
    iszero(dtT) && return st        # dt=0 is a true no-op: advance nothing, and (critically) do
    #                                 NOT consume the step==0 one-time leapfrog-priming guard.
    return _step_hybrid!(st, ps, dtT, NB, nothing)
end

# A velocity source between the field advance and the integer-time moments.
# Plain stepping passes `nothing`; collisional_hybrid.jl owns the only source.
@inline _apply_hybrid_source!(::Nothing, ps, dtT) = nothing

# Body of `step!` for a validated nonzero `dtT`. A source sees positions and B at
# n+1 and velocities at n+1/2; the carried E^{n+1} is then formed from its result.
function _step_hybrid!(
    st::HybridStepper{D,T},
    ps::_HybridPopulations{D,T},
    dtT::T,
    NB::Integer,
    source,
) where {D,T}
    next_time, next_step = _validated_step_clock(st.time[], st.step[], dtT)
    count = _validate_hybrid_species(ps)
    _validate_pressure_step_state(st, ps)
    _save_pressure_transaction!(st, ps)
    completed = false
    try
        _resize_hybrid_particle_workspaces!(st, count)
        g = st.g
        h = dtT / 2
        kick_dt = st.step[] == 0 ? dtT : st.last_dt[] / 2 + h
        # 1–2. Push populations and combine all frozen n+1/2 charge moments.
        _push_hybrid_midpoint_moments!(st, ps, dtT, kick_dt)
        _advance_hybrid_fields!(st, dtT, NB)
        _apply_hybrid_source!(source, ps, dtT)

        # 4. recompute carried E = E^{n+1} (n^{n+1}, u_i^{n+1/2}, B^{n+1}), then re-center u_i to n+1
        #    via a predictor half-kick so the carried E is 2nd-order accurate (see _recenter_carried_E!).
        _hybrid_integer_moments!(st, ps)
        _carried_ohms_law!(st)
        _recenter_carried_E!(st, ps, dtT)
        _validate_pressure_electric!(st)

        st.time[] = next_time
        st.step[] = next_step
        st.last_dt[] = dtT
        completed = true
        return st
    finally
        completed || _restore_pressure_transaction!(st, ps)
    end
end

# ---------------------------------------------------------------- diagnostics

@inline _kinetic_factor_certified(value, lower, upper) =
    iszero(value) | ((lower <= abs(value)) & (abs(value) <= upper))

function _kinetic_energy_bounds(::Type{T}, count::Int128) where {T}
    count_bits = count <= 1 ? 0 : 128 - leading_zeros(count - 1)
    upper = fld(exponent(floatmax(T)) - count_bits - 4, 4)
    positive = min(upper, fld(-exponent(floatmin(T)) - 1, 4))
    mixed = min(upper, fld(-exponent(floatmin(T)) - 4precision(T) + 3, 4))
    positive_bounds =
        positive < 0 ? (T(Inf), zero(T)) : (ldexp(one(T), -positive), ldexp(one(T), positive))
    mixed_bounds = mixed < 0 ? (T(Inf), zero(T)) : (ldexp(one(T), -mixed), ldexp(one(T), mixed))
    return (positive_bounds..., mixed_bounds...)
end

# One ordinary reduction owns scalar and collection arithmetic. Certificates
# exclude intermediate range loss; uncertain finite cases use one complete sum.
function _kinetic_energy_probe(ps::ParticleSet{D,T}, bounds, ::Val{check}) where {D,T,check}
    positive_ok = mixed_ok = true
    has_positive = has_negative = false
    vx, vy, vz = ps.v
    if check
        lo, hi, mixed_lo, mixed_hi = bounds
        positive_ok = _kinetic_factor_certified(ps.m, lo, hi)
        mixed_ok = _kinetic_factor_certified(ps.m, mixed_lo, mixed_hi)
        nonzero_mass, negative_mass = !iszero(ps.m), signbit(ps.m)
        # Boolean reductions may vectorize; the floating reduction below keeps
        # its original order. All inputs are read-only, including shared views.
        @inbounds @simd for p in eachindex(ps.weight)
            w, x, y, z = ps.weight[p], vx[p], vy[p], vz[p]
            positive_ok &=
                _kinetic_factor_certified(w, lo, hi) &
                _kinetic_factor_certified(x, lo, hi) &
                _kinetic_factor_certified(y, lo, hi) &
                _kinetic_factor_certified(z, lo, hi)
            mixed_ok &=
                _kinetic_factor_certified(w, mixed_lo, mixed_hi) &
                _kinetic_factor_certified(x, mixed_lo, mixed_hi) &
                _kinetic_factor_certified(y, mixed_lo, mixed_hi) &
                _kinetic_factor_certified(z, mixed_lo, mixed_hi)
            contributes = nonzero_mass & !iszero(w) & !(iszero(x) & iszero(y) & iszero(z))
            same_sign = negative_mass == signbit(w)
            has_positive |= contributes & same_sign
            has_negative |= contributes & !same_sign
        end
    end
    s = zero(T)
    @inbounds for p in eachindex(ps.weight)
        s += ps.weight[p] * (vx[p]^2 + vy[p]^2 + vz[p]^2)
    end
    return (T(0.5) * ps.m * s, positive_ok, mixed_ok, has_positive, has_negative)
end

@inline function _kinetic_normal_product(a::T, b::T, product::T) where {T}
    # A rounded value at floatmin can have a subnormal exact product. The
    # factor of two certifies the exact product's normal range as well.
    return isfinite(product) && (iszero(a) || iszero(b) || abs(product) >= T(2) * floatmin(T))
end

function _kinetic_energy_absolute_probe(species, ::Type{T}, ::Val{check}) where {T,check}
    absolute = zero(T)
    normal = true
    for ps in species
        sum_abs = zero(T)
        vx, vy, vz = ps.v
        @inbounds for p in eachindex(ps.weight)
            x, y, z = vx[p], vy[p], vz[p]
            x2, y2, z2 = x^2, y^2, z^2
            speed2 = x2 + y2 + z2
            weight = abs(ps.weight[p])
            weighted = weight * speed2
            if check
                normal &=
                    _kinetic_normal_product(x, x, x2) &&
                    _kinetic_normal_product(y, y, y2) &&
                    _kinetic_normal_product(z, z, z2) &&
                    isfinite(speed2) &&
                    _kinetic_normal_product(weight, speed2, weighted)
            end
            sum_abs += weighted
        end
        mass = abs(ps.m)
        half_mass = T(0.5) * mass
        contribution = half_mass * sum_abs
        if check
            normal &=
                isfinite(sum_abs) &&
                _kinetic_normal_product(T(0.5), mass, half_mass) &&
                _kinetic_normal_product(half_mass, sum_abs, contribution)
        end
        absolute += contribution
    end
    return absolute, normal && isfinite(absolute)
end

_kinetic_energy_absolute(species, ::Type{T}) where {T} =
    first(_kinetic_energy_absolute_probe(species, T, Val(false)))

# Wide-certified factors keep products normal. Subnormal additions are exact
# on the minimum-subnormal lattice. A contribution has at most count+4 rounded
# operations; the positive absolute evaluation has the same bound.
function _kinetic_energy_error_bound(absolute::T, count::Int128) where {T}
    k = count + 4
    k <= Int128(1) << (precision(T) - 2) || return T(Inf)
    delta = ldexp(T(k), -precision(T))
    factor = nextfloat(delta / (one(T) - T(2) * delta))
    return nextfloat(factor * absolute)
end

function _kinetic_energy_normal_certified(species, value::T, count::Int128) where {T}
    count + 4 <= Int128(1) << (precision(T) - 2) || return false
    absolute = _kinetic_energy_absolute(species, T)
    error = _kinetic_energy_error_bound(absolute, count)
    return prevfloat(abs(value) - error) >= floatmin(T)
end

# Direct product checks cover ordinary inputs outside the uniform factor box.
# Signed population products can underflow through cancellation. Their combined
# propagated rounding error is at most 2*count*qmin/3 under the same count limit;
# count*qmin is exactly representable and gives a conservative extra allowance.
function _kinetic_energy_operations_certified(species, value::T, count::Int128) where {T}
    count + 4 <= Int128(1) << (precision(T) - 2) || return false
    absolute, normal = _kinetic_energy_absolute_probe(species, T, Val(true))
    normal || return false
    error = nextfloat(_kinetic_energy_error_bound(absolute, count) + T(count) * nextfloat(zero(T)))
    return prevfloat(abs(value) - error) >= floatmin(T) && isfinite(nextfloat(abs(value) + error))
end

function _kinetic_energy_exact(species, ::Type{T}) where {T}
    total = zero(Rational{BigInt})
    for ps in species
        mass = _exact_rational(ps.m)
        iszero(mass) && continue
        for p in eachindex(ps.weight)
            weight = _exact_rational(ps.weight[p])
            iszero(weight) && continue
            speed2 = sum(_exact_rational(ps.v[c][p])^2 for c = 1:3)
            total += mass * weight * speed2 / 2
        end
    end
    return _from_exact(T, total)
end

function _kinetic_energy_checked(species, ::Type{T}, ::Val{single}) where {T,single}
    count = Int128(0)
    for ps in species
        # Both collection and population lengths fit Int; their complete
        # product fits Int128, including the largest 64-bit Int lengths.
        count += length(ps.weight)
    end
    bounds = _kinetic_energy_bounds(T, count)
    energy = zero(T)
    positive_ok = mixed_ok = true
    has_positive = has_negative = false
    for ps in species
        value, positive, mixed, plus, minus = _kinetic_energy_probe(ps, bounds, Val(true))
        energy = single ? value : energy + value
        positive_ok &= positive
        mixed_ok &= mixed
        has_positive |= plus
        has_negative |= minus
    end
    # Empty or identically zero contributions keep the ordinary signed zero.
    (iszero(count) || (!(has_positive || has_negative) && iszero(energy))) && return energy
    (mixed_ok || (positive_ok && !(has_positive && has_negative))) && return energy
    positive_ok && _kinetic_energy_normal_certified(species, energy, count) && return energy
    !positive_ok && _kinetic_energy_operations_certified(species, energy, count) && return energy
    finite = all(species) do ps
        isfinite(ps.m) && all(isfinite, ps.weight) && all(a -> all(isfinite, a), ps.v)
    end
    # Rational conversion is eligible only when every finite check is known true.
    return finite === true ? _kinetic_energy_exact(species, T) : energy
end

"""
    kinetic_energy(ps) -> Real

Total ion kinetic energy `Σ_p ½ m w |v|²`. For Float16/32/64, range-sensitive
finite factors are retained through the complete expression before final
rounding. Inputs certified free of range loss keep their existing summation
order. Signed weights and masses are supported; nonfinite inputs retain their
arithmetic behavior.
"""
function kinetic_energy(ps::ParticleSet{D,T}) where {D,T}
    if T <: Union{Float16,Float32,Float64}
        return _kinetic_energy_checked((ps,), T, Val(true))
    end
    return first(_kinetic_energy_probe(ps, nothing, Val(false)))
end

"""
    kinetic_energy(species::AbstractVector) -> Real

Total kinetic energy of all populations; a typed empty collection contributes
zero. Range-sensitive finite Float16/32/64 inputs retain the whole collection
before rounding, including cancellation between overflowing population energies.
"""
function kinetic_energy(species::AbstractVector{<:ParticleSet{D,T}}) where {D,T}
    if T <: Union{Float16,Float32,Float64}
        return _kinetic_energy_checked(species, T, Val(false))
    end
    energy = zero(T)
    for ps in species
        energy += kinetic_energy(ps)
    end
    return energy
end

# Cold energy arithmetic retains all finite particle factors through final rounding.

function _relativistic_rational_sqrt(value::Rational{BigInt})
    n, d = numerator(value), denominator(value)
    sn, sd = isqrt(n), isqrt(d)
    return sn * sn == n && sd * sd == d ? sn // sd : nothing
end

function _combine_relativistic_radical!(
    radicands::Vector{Rational{BigInt}},
    coefficients::Vector{Rational{BigInt}},
    r::Rational{BigInt},
    a::Rational{BigInt},
)
    iszero(a) && return zero(Rational{BigInt})
    root = _relativistic_rational_sqrt(r)
    root === nothing || return a * root
    for i in eachindex(radicands)
        ratio_root = _relativistic_rational_sqrt(r / radicands[i])
        if ratio_root !== nothing
            coefficients[i] += a * ratio_root
            return zero(Rational{BigInt})
        end
    end
    push!(radicands, r)
    push!(coefficients, a)
    return zero(Rational{BigInt})
end

function _finish_relativistic_energy(
    ::Type{T},
    rational::Rational{BigInt},
    radicands,
    coefficients,
) where {T}
    active = findall(!iszero, coefficients)
    isempty(active) && return _from_exact(T, rational)
    bits = precision(T) + 16
    canonical = length(active) < 2
    while true
        scale = BigInt(1) << bits
        lower, upper = rational, rational
        for i in active
            r, a = radicands[i], coefficients[i]
            # floor(sqrt(r)*2^bits), obtained without a floating-point operation.
            root_floor = isqrt(div(numerator(r) << (2bits), denominator(r)))
            lo, hi = root_floor // scale, (root_floor + 1) // scale
            if a > 0
                lower += a * lo
                upper += a * hi
            else
                lower += a * hi
                upper += a * lo
            end
        end
        rounded_lower = _from_exact(T, lower)
        rounded_upper = _from_exact(T, upper)
        if isequal(rounded_lower, rounded_upper)
            return rounded_lower
        end
        if !canonical && (
            lower <= 0 <= upper ||
            rounded_lower == rounded_upper ||
            (isfinite(rounded_lower) && nextfloat(rounded_lower) == rounded_upper)
        )
            # Only exact zero or an exact rational midpoint can prevent eventual
            # agreement. Resolve square-class relations when the bounds have
            # narrowed to one rounding boundary or enclose zero, not for every
            # ordinary sum. BigFloat zero needs this explicit symbolic path.
            combined_r, combined_a = Rational{BigInt}[], Rational{BigInt}[]
            for i in active
                rational += _combine_relativistic_radical!(
                    combined_r,
                    combined_a,
                    radicands[i],
                    coefficients[i],
                )
            end
            radicands, coefficients = combined_r, combined_a
            active = findall(!iszero, coefficients)
            isempty(active) && return _from_exact(T, rational)
            canonical = true
        end
        bits *= 2
    end
end

function _kinetic_energy_relativistic_exact(ps::ParticleSet{D,T}, c::Real) where {D,T}
    cT = _require_finite_positive_real("c", c, T)
    c2 = _exact_rational(cT)^2
    mc2 = _exact_rational(ps.m) * c2
    rational = zero(Rational{BigInt})
    radicands, coefficients = Rational{BigInt}[], Rational{BigInt}[]
    indices = Dict{Rational{BigInt},Int}()
    for p in eachindex(ps.weight)
        all(isfinite(ps.v[d][p]) for d = 1:3) ||
            throw(ArgumentError("particle $p velocity must be finite"))
        v2 = sum(_exact_rational(ps.v[d][p])^2 for d = 1:3)
        v2 < c2 || throw(ArgumentError("particle $p has |v| ≥ c"))
        iszero(v2) && continue
        a = mc2 * _exact_rational(ps.weight[p])
        iszero(a) && continue
        rational -= a
        r = c2 / (c2 - v2)
        index = get(indices, r, 0)
        if index != 0
            coefficients[index] += a
        else
            root = _relativistic_rational_sqrt(r)
            if root === nothing
                push!(radicands, r)
                push!(coefficients, a)
                indices[r] = length(radicands)
            else
                rational += a * root
            end
        end
    end
    return _finish_relativistic_energy(T, rational, radicands, coefficients)
end

@inline _relativistic_range_lost(value::T) where {T<:AbstractFloat} =
    !isfinite(value) | (!iszero(value) & (abs(value) < floatmin(T)))
@inline _relativistic_range_lost(value::Union{Float16,Float32,Float64}) =
    !isfinite(value) | issubnormal(value)
@inline _relativistic_product_lost(a, b, value) =
    _relativistic_range_lost(value) | (iszero(value) & !iszero(a) & !iszero(b))

# Nonfinite weights keep their original multiplication-before-division order.
# Only the uncertain norm and root are resolved from the exact stored components.
function _relativistic_beta_root(vx, vy, vz, c, particle)
    (isfinite(vx) && isfinite(vy) && isfinite(vz)) ||
        throw(ArgumentError("particle $particle velocity must be finite"))
    c2 = _exact_rational(c)^2
    v2 = _exact_rational(vx)^2 + _exact_rational(vy)^2 + _exact_rational(vz)^2
    v2 < c2 || throw(ArgumentError("particle $particle has |v| ≥ c"))
    ratio = (c2 - v2) / c2
    exact_root = _relativistic_rational_sqrt(ratio)
    T = typeof(c)
    root =
        exact_root === nothing ?
        _finish_relativistic_energy(T, zero(Rational{BigInt}), [ratio], [one(Rational{BigInt})]) :
        _from_exact(T, exact_root)
    return _from_exact(T, v2 / c2), root
end

@inline _relativistic_magnitude_certified(value, lower, upper) =
    iszero(value) | ((lower <= abs(value)) & (abs(value) <= upper))

function _relativistic_range_certificate(ps::ParticleSet{D,T}, c::T) where {D,T}
    T <: Union{Float32,Float64} || return false
    p = precision(T)
    nbits = ndigits(max(1, length(ps.weight)); base = 2)
    # Eight bits of slack bound the contribution lattice, partial sums and
    # final mass/c factor away from both normal-range limits.
    k = min(fld(-exponent(floatmin(T)) - p - 8, 8), fld(exponent(floatmax(T)) - p - nbits - 8, 4))
    k > 0 || return false
    lower, upper = ldexp(one(T), -k), ldexp(one(T), k)
    (
        _relativistic_magnitude_certified(c, lower, upper) &&
        _relativistic_magnitude_certified(ps.m, lower, upper)
    ) || return false
    vx, vy, vz = ps.v
    certified = true
    @inbounds @simd for i in eachindex(ps.weight)
        certified &=
            _relativistic_magnitude_certified(vx[i], lower, upper) &
            _relativistic_magnitude_certified(vy[i], lower, upper) &
            _relativistic_magnitude_certified(vz[i], lower, upper) &
            _relativistic_magnitude_certified(ps.weight[i], lower, upper)
    end
    return certified
end

"""
    kinetic_energy_relativistic(ps, c) -> Real

Total relativistic kinetic energy `Σ_p w m (γ_p − 1) c²` with
`γ_p = 1/√(1 − |v_p|²/c²)`, for `relativistic = true` EMPIC runs (`ParticleSet`
stores the velocity `v`; the relativistic pusher works in `u = γv` internally).
In the Newtonian limit `|v| ≪ c` this reduces to [`kinetic_energy`](@ref)`(ps)`
— the per-particle difference is `O(w m |v|⁴/c²)`, i.e. a relative `O(v²/c²)`.
Throws `ArgumentError` if any particle has `|v| ≥ c`: the relativistic pusher
keeps `|v| < c`, so such a particle is corrupted state, not a valid energy.
Uncertain speed comparisons use the exact stored components. For finite mass
and weights, intermediate overflow or underflow triggers evaluation of the
complete weighted expression before final rounding. Ordinary inputs retain
their original arithmetic order; signed and nonfinite factors remain visible.
"""
function kinetic_energy_relativistic(ps::ParticleSet{D,T}, c::Real) where {D,T}
    cT = _require_finite_positive_real("c", c, T)
    if _relativistic_range_certificate(ps, cT)
        return _kinetic_energy_relativistic(ps, cT, Val(false))
    end
    return _kinetic_energy_relativistic(ps, cT, Val(true))
end

function _kinetic_energy_relativistic(
    ps::ParticleSet{D,T},
    cT::T,
    ::Val{check_range},
) where {D,T,check_range}
    c2 = cT * cT
    finite_factors = !check_range || (isfinite(ps.m) && all(isfinite, ps.weight))
    if check_range && finite_factors && (_relativistic_range_lost(c2) || iszero(c2))
        return _kinetic_energy_relativistic_exact(ps, cT)
    end
    s = zero(T)
    range_lost = false
    vx, vy, vz = ps.v
    @inbounds for p in eachindex(ps.weight)
        x, y, z = vx[p], vy[p], vz[p]
        x2, y2, z2 = x * x, y * y, z * z
        v2 = x2 + y2 + z2
        # Resolve either side of the uncertain boundary. Rounded velocity
        # squares can falsely accept a superluminal stored tuple as well as
        # falsely reject a valid one. The pusher has a separate, narrower range.
        uncertain =
            _relativistic_range_lost(c2) ||
            iszero(c2) ||
            !isfinite(v2) ||
            v2 >= c2 ||
            c2 - v2 <= T(8) * eps(T) * c2
        if uncertain
            finite_factors && return _kinetic_energy_relativistic_exact(ps, cT)
            beta2, r = _relativistic_beta_root(x, y, z, cT, p)
        else
            beta2 = v2 / c2
            r = sqrt(one(T) - beta2)
        end
        if check_range
            range_lost |=
                _relativistic_product_lost(x, x, x2) |
                _relativistic_product_lost(y, y, y2) |
                _relativistic_product_lost(z, z, z2) |
                _relativistic_range_lost(beta2) |
                (iszero(beta2) & !iszero(v2))
        end
        weighted = ps.weight[p] * beta2
        denominator = r * (one(T) + r)
        contribution = weighted / denominator
        next_sum = s + contribution
        if check_range
            range_lost |=
                _relativistic_product_lost(ps.weight[p], beta2, weighted) |
                _relativistic_range_lost(contribution) |
                (iszero(contribution) & !iszero(weighted)) |
                _relativistic_range_lost(next_sum)
        end
        s = next_sum
    end
    mc2 = ps.m * c2
    value = mc2 * s
    if check_range &&
       finite_factors &&
       (
           range_lost ||
           _relativistic_product_lost(ps.m, c2, mc2) ||
           _relativistic_product_lost(mc2, s, value)
       )
        return _kinetic_energy_relativistic_exact(ps, cT)
    end
    return value
end

"""
    magnetic_energy(B, g) -> Real

Magnetic energy `∫ ½|B|² dV` in normalized units. Finite Float32 fields use
certified wide accumulation and exact completion before rounding, sharing the
arithmetic owner of [`em_field_energy`](@ref). Other precisions retain their
ordinary summation order and range repair.
"""
function magnetic_energy(B::NTuple{3,<:Array{T,D}}, g::FourierGrid{D,T}) where {D,T}
    _require_grid_tuple(:B, B, g)
    if T === Float32
        value = _field_energy_float32((), B, 1.0, g.dx)
        value === nothing || return value
    end
    s = zero(T)
    minimum_b2 = T(Inf)
    @inbounds for c = 1:3, I in eachindex(B[c])
        value = B[c][I]
        square = value^2
        minimum_b2 = min(minimum_b2, ifelse(iszero(value), T(Inf), square))
        s += square
    end
    return _em_field_energy_finish((), B, one(T), g.dx, s, T(Inf), minimum_b2)
end

# For nonnegative inputs, a finite rounded sum is at least every input, hence
# the exact sum is at most count*pressure_sum. Numerator factors are bounded
# above by powers of two. The predecessor of rounded gamma-1 bounds its exact
# positive value from below. Int128 retains the sum of machine-sized exponents.
function _pressure_energy_upper_exponent(pressure_sum::T, γ::T, dx, count::Int) where {T}
    count_exponent = 8sizeof(Int) - leading_zeros(count - 1)
    upper = Int128(exponent(pressure_sum)) + 1 + count_exponent
    for spacing in dx
        upper += Int128(exponent(spacing)) + 1
    end
    return upper - Int128(exponent(prevfloat(γ - one(T))))
end

# Pressure and gamma have been validated by the public diagnostic. Preserve
# ordinary summation/order; recover the complete integral on intermediate range
# loss or when the positive upper bound cannot exclude final overflow.
function _stored_pressure_internal_energy(pe::AbstractArray{T}, γ::T, dx::NTuple{D,T}) where {D,T}
    pressure_sum = sum(pe)
    denominator = γ - one(T)
    unscaled = pressure_sum / denominator
    volume = prod(dx)
    result = unscaled * volume
    lost =
        !isfinite(pressure_sum) ||
        _spectral_quotient_lost_range(pressure_sum, denominator, unscaled) ||
        _spectral_product_lost_range(unscaled, volume, result)
    # Contiguous products cover either multiplication association in D <= 3.
    for first = 1:length(dx)
        partial = one(T)
        for last = first:length(dx)
            product = partial * dx[last]
            lost |= _spectral_product_lost_range(partial, dx[last], product)
            partial = product
        end
    end
    if !lost && !iszero(pressure_sum)
        upper = _pressure_energy_upper_exponent(pressure_sum, γ, dx, length(pe))
        lost = upper > exponent(floatmax(T))
    end
    if lost
        total = zero(Rational{BigInt})
        for pressure in pe
            total += _exact_rational(pressure)
        end
        integral = total * prod(_exact_rational, dx) / (_exact_rational(γ) - 1)
        return _from_exact(T, integral)
    end
    return result
end

"""
    electron_internal_energy(n, closure, g; pe=nothing) -> Real
    electron_internal_energy(n, B, closure, g; pe=nothing) -> Real

Electron closure energy integrated over the grid:

  - evolving scalar pressure: `∫ pe/(γ−1) dV`, requiring the stored state through
    `pe=stepper.fields.pe`; density alone cannot reconstruct this reservoir;
  - polytropic `γ ≠ 1`: `∫ p_e/(γ−1) dV` from the algebraic closure;
  - `γ = 1` (isothermal): electron free energy `T_e ∫ n ln n dV`, with cells
    at `n ≤ 0` contributing zero;
  - CGL (gyrotropic): `∫ (p_⊥ + p_∥/2) dV`, requiring the B-taking method.
    The density-only method retains its existing `NaN` result for CGL.

For finite Float32/64 isothermal data and temperature, states whose positive
densities are all within `1/256` of one use bounded complete free-energy
evaluation. Other finite states use a compensated Float64 reduction. Unresolved
bounds or range-sensitive operations use the complete stored-input expression
before final rounding. The compensated path does not guarantee single rounding.
Nonpositive densities retain their existing zero contribution.

For finite Float32/64 algebraic polytropic inputs, branch selection retains the
stored closure gamma, including when closure and grid precisions differ.
Ordinary `gamma != 1` states use compensated accumulation; range-sensitive
cases retain the complete density ratio, power, coefficient and cell spacings
before final rounding. This ordinary path does not guarantee single rounding.
Exact stored `gamma == 1` uses the complete logarithmic expression with the
coefficient `pe0/n0`, without rounding that ratio first.

For finite Float32/64 CGL inputs and invariant coefficients, range-sensitive
cases retain the stored density, field components, coefficients and individual
cell spacings through the complete integral before final rounding. Inputs
certified free of range loss retain the existing floating arithmetic. At zero
field with nonzero parallel pressure, the existing singular `Inf`/`NaN` result
is preserved. This diagnostic uses the supplied field magnitude, without the
floor used by the pressure-force calculation.

Algebraic closures reject an alternate `pe` argument. This function integrates
supplied state without evolving it. See [`energy_budget`](@ref) for continuum
model qualifications and particle/field time levels; this diagnostic does not
imply exact discrete energy conservation.
"""
function electron_internal_energy(
    n::Array{T,D},
    closure::ElectronClosure,
    g::FourierGrid{D,T};
    pe = nothing,
) where {D,T}
    _require_grid_array(:n, n, g)
    if _evolves_pressure(closure)
        pe === nothing && throw(
            ArgumentError(
                "$(nameof(typeof(closure))) energy requires explicit stored pressure via pe=",
            ),
        )
        pe isa AbstractArray{T,D} ||
            throw(ArgumentError("stored pressure must match grid precision and dimension"))
        _validate_evolving_pressure(pe, g)
        γ = _require_finite_real("γ", closure.γ, T)
        γ > one(T) ||
            throw(ArgumentError("evolving electron γ must be greater than one at grid precision"))
        return _stored_pressure_internal_energy(pe, γ, g.dx)
    end
    pe === nothing ||
        throw(ArgumentError("algebraic electron energy does not accept an alternate pe"))
    # a gyrotropic (CGL) closure needs |B| for its internal energy — use the
    # (n, B, closure, g) method; degrade to NaN here rather than throwing downstream.
    is_anisotropic(closure) && return T(NaN)
    if T <: Union{Float32,Float64} && closure isa IsothermalElectrons{<:Union{Float32,Float64}}
        return _isothermal_energy(n, closure.Te, g.dx)
    end
    if T <: Union{Float32,Float64} && closure isa PolytropicElectrons{<:Union{Float32,Float64}}
        if all(isfinite, (closure.pe0, closure.n0, closure.γ)) &&
           closure.pe0 >= 0 &&
           closure.n0 > 0 &&
           closure.γ > 0 &&
           all(isfinite, n) &&
           all(isfinite, g.dx)
            return _polytropic_energy(n, closure, g.dx)
        end
    end
    γ = T(closure_gamma(closure))
    pe = similar(n)
    electron_pressure!(pe, n, closure)
    if γ == one(T)
        # isothermal free energy F_e = T_e ∫ n ln n dV = ∫ p_e ln n dV (p_e = T_e n)
        s = zero(T)
        @inbounds for I in eachindex(pe, n)
            n[I] > zero(T) && (s += pe[I] * log(n[I]))
        end
        return s * prod(g.dx)
    end
    return sum(pe) / (γ - one(T)) * prod(g.dx)
end

# scalar closures do not need B — accept it for interface symmetry with the CGL method
function electron_internal_energy(
    n::Array{T,D},
    B::NTuple{3,<:Array{T,D}},
    closure::ElectronClosure,
    g::FourierGrid{D,T};
    pe = nothing,
) where {D,T}
    _require_grid_tuple(:B, B, g)
    return electron_internal_energy(n, closure, g; pe)
end

# The ordinary CGL reduction and its error certificate share one traversal.
@inline function _cgl_energy_normal_product(a, b, product::T) where {T}
    return isfinite(product) && (iszero(a) || iszero(b) || abs(product) >= T(2) * floatmin(T))
end
@inline function _cgl_energy_normal_conversion(a, converted::T) where {T}
    return isfinite(converted) && (iszero(a) || abs(converted) >= T(2) * floatmin(T))
end
@inline function _cgl_energy_normal_quotient(a, b, quotient::T) where {T}
    return isfinite(b) &&
           !iszero(b) &&
           isfinite(quotient) &&
           (iszero(a) || abs(quotient) >= T(2) * floatmin(T))
end

function _cgl_energy_probe(n::Array{T,D}, B, closure, dx, ::Val{check}) where {T,D,check}
    s = zero(T)
    absolute = largest = zero(T)
    normal = true
    active = false
    nonnegative = !check || (closure.Cperp >= 0 && closure.Cpar >= 0)
    @inbounds for I in eachindex(n, B[1], B[2], B[3])
        bx, by, bz = B[1][I], B[2][I], B[3][I]
        x2, y2, z2 = bx^2, by^2, bz^2
        norm2 = x2 + y2 + z2
        Bmag = sqrt(norm2)
        perpendicular = T(cgl_pperp(closure, n[I], Bmag))
        parallel = T(cgl_ppar(closure, n[I], Bmag))
        half_parallel = T(0.5) * parallel
        s += perpendicular + half_parallel
        if check
            absolute += abs(perpendicular) + abs(half_parallel)
            largest = max(largest, abs(perpendicular), abs(half_parallel))
            nonnegative &= n[I] >= zero(T)
            zero_cell =
                iszero(n[I]) ||
                (iszero(closure.Cperp) && iszero(closure.Cpar)) ||
                (iszero(bx) && iszero(by) && iszero(bz) && iszero(closure.Cpar))
            zero_cell && continue
            active = true
            C = typeof(closure.Cperp)
            density, field = C(n[I]), C(Bmag)
            normal &=
                _cgl_energy_normal_product(bx, bx, x2) &&
                _cgl_energy_normal_product(by, by, y2) &&
                _cgl_energy_normal_product(bz, bz, z2) &&
                isfinite(norm2) &&
                _cgl_energy_normal_conversion(n[I], density) &&
                _cgl_energy_normal_conversion(Bmag, field)
            if !iszero(closure.Cperp)
                perp_product = closure.Cperp * density
                raw_perpendicular = perp_product * field
                normal &=
                    _cgl_energy_normal_product(closure.Cperp, density, perp_product) &&
                    _cgl_energy_normal_product(perp_product, field, raw_perpendicular) &&
                    _cgl_energy_normal_conversion(raw_perpendicular, perpendicular)
            end
            if !iszero(closure.Cpar)
                par_product = closure.Cpar * density
                par_square = par_product * density
                par_cube = par_square * density
                denominator = field * field
                raw_parallel = par_cube / denominator
                normal &=
                    _cgl_energy_normal_product(closure.Cpar, density, par_product) &&
                    _cgl_energy_normal_product(par_product, density, par_square) &&
                    _cgl_energy_normal_product(par_square, density, par_cube) &&
                    _cgl_energy_normal_product(field, field, denominator) &&
                    _cgl_energy_normal_quotient(par_cube, denominator, raw_parallel) &&
                    _cgl_energy_normal_conversion(raw_parallel, parallel) &&
                    _cgl_energy_normal_product(T(0.5), parallel, half_parallel)
            end
        end
    end
    volume = prod(dx)
    result = s * volume
    if check
        unscaled_absolute = absolute
        absolute *= abs(volume)
        normal &=
            isfinite(s) &&
            isfinite(unscaled_absolute) &&
            _cgl_energy_normal_product(s, volume, result) &&
            _cgl_energy_normal_product(unscaled_absolute, abs(volume), absolute)
        for first in eachindex(dx)
            partial = one(T)
            for last = first:length(dx)
                product = partial * dx[last]
                normal &= _cgl_energy_normal_product(partial, dx[last], product)
                partial = product
            end
        end
    end
    return result, absolute, normal, active, largest, nonnegative
end

# Normal operations bound each component's exact/rounded ratio in both
# directions by 2^19. At most two nonnegative components contribute per cell.
# Individual spacing exponents retain the complete volume in both bounds.
function _cgl_energy_positive_certified(largest::T, dx, count::Int) where {T}
    (count > 0 && isfinite(largest) && largest > zero(T)) || return false
    count_exponent = 128 - leading_zeros(Int128(count) - 1)
    lower = Int128(exponent(largest)) - 19
    upper = Int128(exponent(largest)) + 21 + count_exponent
    for spacing in dx
        (isfinite(spacing) && !iszero(spacing)) || return false
        e = Int128(exponent(abs(spacing)))
        lower += e
        upper += e + 1
    end
    return lower >= exponent(floatmin(T)) && upper <= exponent(floatmax(T))
end

function _cgl_energy_normal_certified(value::T, absolute::T, count::Int, ::Type{C}) where {T,C}
    p = min(precision(T), precision(C))
    k = Int128(count) + 22
    k <= Int128(1) << (p - 2) || return false
    # A parallel term has at most 19 normal rounding factors, including
    # closure/grid conversions and the squared norm's inverse factors. Add
    # one cell addition, count-1 grid additions, and <=3 volume products.
    # The positive absolute evaluation has the same bound. Signed subnormal
    # additions are exact on the field-precision lattice.
    delta = ldexp(T(k), -p)
    factor = nextfloat(delta / (one(T) - T(2) * delta))
    error = nextfloat(factor * absolute)
    return prevfloat(abs(value) - error) >= floatmin(T) && isfinite(nextfloat(abs(value) + error))
end

function _cgl_energy_exact(n::Array{T,D}, B, closure, dx, ordinary) where {T,D}
    Q = Rational{BigInt}
    perpendicular, parallel = _exact_rational(closure.Cperp), _exact_rational(closure.Cpar)
    volume = prod(_exact_rational, dx)
    rational = zero(Q)
    radicands, coefficients = Q[], Q[]
    indices = Dict{Q,Int}()
    for I in eachindex(n, B[1], B[2], B[3])
        density = _exact_rational(n[I])
        iszero(density) && continue
        norm2 = sum(_exact_rational(component[I])^2 for component in B)
        if iszero(norm2)
            # The existing CGL law is singular for nonzero parallel pressure
            # at zero field; preserve its visible Inf/NaN rather than invent a floor.
            iszero(parallel) || return ordinary
            continue
        end
        rational += parallel * density^3 * volume / (2norm2)
        a = perpendicular * density * volume
        iszero(a) && continue
        index = get(indices, norm2, 0)
        if index != 0
            coefficients[index] += a
        else
            root = _relativistic_rational_sqrt(norm2)
            if root === nothing
                push!(radicands, norm2)
                push!(coefficients, a)
                indices[norm2] = length(radicands)
            else
                rational += a * root
            end
        end
    end
    return _finish_relativistic_energy(T, rational, radicands, coefficients)
end

function _cgl_energy_checked(n::Array{T,D}, B, closure::CGLElectrons{C}, dx) where {T,D,C}
    value, absolute, normal, active, largest, nonnegative =
        _cgl_energy_probe(n, B, closure, dx, Val(true))
    !active && iszero(value) && return value
    if normal
        nonnegative &&
            !iszero(value) &&
            _cgl_energy_positive_certified(largest, dx, length(n)) &&
            return value
        _cgl_energy_normal_certified(value, absolute, length(n), C) && return value
    end
    finite =
        isfinite(closure.Cperp) &&
        isfinite(closure.Cpar) &&
        all(isfinite, n) &&
        all(isfinite, dx) &&
        all(a -> all(isfinite, a), B)
    return finite === true ? _cgl_energy_exact(n, B, closure, dx, value) : value
end

function electron_internal_energy(
    n::Array{T,D},
    B::NTuple{3,<:Array{T,D}},
    closure::CGLElectrons,
    g::FourierGrid{D,T};
    pe = nothing,
) where {D,T}
    pe === nothing ||
        throw(ArgumentError("algebraic electron energy does not accept an alternate pe"))
    _require_grid_array(:n, n, g)
    _require_grid_tuple(:B, B, g)
    # gyrotropic internal energy ½ tr(P_e) = p_⊥ + p_∥/2 with the CGL-slaved
    # pressures p_⊥ = C_⊥ n|B|, p_∥ = C_∥ n³/|B|² (cgl_pperp / cgl_ppar).
    if T <: Union{Float32,Float64} && closure isa CGLElectrons{<:Union{Float32,Float64}}
        return _cgl_energy_checked(n, B, closure, g.dx)
    end
    return first(_cgl_energy_probe(n, B, closure, g.dx, Val(false)))
end

"""
    mode_amplitude(field, g, m::NTuple{D,Int})

Complex Fourier amplitude Σ_g field_g · exp(−i k·x_g) of integer mode `m`
(k_d = 2π m_d / L_d), for tracking a single wave mode over time.
"""
function mode_amplitude(field::Array{T,D}, g::FourierGrid{D,T}, m::NTuple{D,Int}) where {D,T}
    _require_grid_array(:field, field, g)
    acc = zero(Complex{T})
    kk = ntuple(d -> T(2π) * m[d] / g.L[d], D)
    @inbounds for I in CartesianIndices(field)
        t = Tuple(I)
        phase = zero(T)
        for d = 1:D
            phase += kk[d] * (t[d] - 1) * g.dx[d]
        end
        acc += field[I] * cis(-phase)
    end
    return acc
end
