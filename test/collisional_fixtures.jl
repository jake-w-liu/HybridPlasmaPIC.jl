# Shared fixtures for the collisional hybrid test files; each file includes this
# once so it can also run on its own.
using Test, Random, Serialization
using HybridPlasmaPIC

function _composition_fixture(
    ::Val{D},
    ::Type{T},
    R;
    ns = 2,
    count = 2 * 4^D,
    logs = T(0.05),
    temperature = T(0.02),
    magnetic = T(0.1),
    shape = CIC(),
) where {D,T}
    g = FourierGrid(ntuple(_ -> 4, D), ntuple(_ -> T(2), D))
    species = [ParticleSet{D,T}(count; q = T(a), m = T(a + 1)) for a = 1:ns]
    for (a, p) in enumerate(species)
        for j = 1:count, d = 1:D
            p.x[d][j] = T(mod(div(j - 1, 4^(d - 1)), 4) + 0.25 + 0.1a) * T(0.5)
        end
        for j = 1:count, c = 1:3
            p.v[c][j] = T(0.1sin(j + a + c))
        end
        fill!(p.weight, T(2^D / max(count, 1)))
    end
    inner = HybridStepper(g, HybridModel(IsothermalElectrons(temperature)), shape, ns * count)
    fill!(inner.fields.B[1], magnetic)
    st = CollisionalHybridStepper(
        inner,
        species;
        dt = T(0.001),
        rng = R(712),
        coulomb_log = logs,
        epsilon0 = one(T),
        weight_to_number = one(T),
    )
    return st, species
end

# Finite inputs can yield an unrepresentable sampled velocity. No injection hook
# or custom RNG is needed to exercise failure after field and random advancement.
function _composition_sampled_failure_fixture(::Type{T}, R, seed) where {T}
    g = FourierGrid((1,), (one(T),))
    charge = floatmax(T) / T(2)
    p = ParticleSet{1,T}(2; q = charge, m = charge)
    fill!(p.x[1], T(0.25))
    fill!(p.weight, T(0.001) / charge)
    fill!(p.v[1], T(0.95) * floatmax(T))
    p.v[2] .= (T(0.1) * floatmax(T), -T(0.1) * floatmax(T))
    inner = HybridStepper(g, HybridModel(IsothermalElectrons(zero(T)); nfloor = 1.0), NGP(), 2)
    st = CollisionalHybridStepper(
        inner,
        [p];
        dt = T(0.001),
        rng = R(seed),
        coulomb_log = one(T),
        epsilon0 = nextfloat(zero(T)),
        weight_to_number = one(T),
        u_floor = one(T),
    )
    return st, [p]
end

function _composition_rng_state(rng)
    return deepcopy(
        NamedTuple{fieldnames(typeof(rng))}(ntuple(i -> getfield(rng, i), fieldcount(typeof(rng)))),
    )
end

function _composition_snapshot(st, species)
    H = HybridPlasmaPIC
    return deepcopy((
        particles = [H._checkpoint_particle_state(p) for p in species],
        fields = H._hybrid_field_state(st.inner.fields),
        floor = st.inner.fields.floor_count[],
        clock = (st.inner.time[], st.inner.step[], st.inner.last_dt[]),
        ledger = (st.source_time, st.source_count, st.phase),
        rng = _composition_rng_state(st.rng),
    ))
end

function _composition_physical(st, species)
    s = _composition_snapshot(st, species)
    return (; s.particles, s.fields, s.floor, s.clock, s.rng)
end

function _composition_component!(st, species)
    H = HybridPlasmaPIC
    step!(st.inner, species, st.dt)
    # Independent reference: the plain stepper, the public collision call, then
    # the canonical moment/Ohm/carried-E recomputation from collided velocities.
    if any(!iszero, st.coulomb_log)
        collide_coulomb!(
            species,
            st.inner.g,
            st.dt;
            coulomb_log = st.coulomb_log,
            epsilon0 = st.epsilon0,
            weight_to_number = st.weight_to_number,
            u_floor = st.u_floor,
            rng = st.rng,
            workspace = st.workspace,
        )
        H._hybrid_integer_moments!(st.inner, species)
        ohms_law!(st.inner.fields, st.inner.model, st.inner.g)
        H._recenter_carried_E!(st.inner, species, st.dt)
    end
end

# Odd grids, NGP, cold electrons and one particle cluster near every node per
# species (a particles of species a), so every particle gathers its own node's E.
function _composition_feedback_fixture(::Val{D}, ::Type{T}, R; ns = 2, strong = false) where {D,T}
    dims = ntuple(d -> d == 1 ? 5 : 3, D)
    g = FourierGrid(dims, ntuple(_ -> T(2pi), D))
    species = [ParticleSet{D,T}(a * prod(dims); q = T(a), m = T(a + 1)) for a = 1:ns]
    for (a, p) in enumerate(species), (node, I) in enumerate(CartesianIndices(dims)), j = 1:a
        index = (node - 1) * a + j
        for d = 1:D
            p.x[d][index] = (T(I[d] - 1) + T(0.2 + 0.01a + 0.005j)) * g.dx[d]
        end
        for c = 1:3
            p.v[c][index] = T(0.1sin(node + a + j + c))
        end
        p.weight[index] = prod(g.dx)
    end
    inner =
        HybridStepper(g, HybridModel(IsothermalElectrons(zero(T))), NGP(), sum(nparticles, species))
    for I in CartesianIndices(dims)
        x = T(I[1] - 1) * g.dx[1]
        y = D >= 2 ? T(I[2] - 1) * g.dx[2] : zero(T)
        z = D >= 3 ? T(I[3] - 1) * g.dx[3] : zero(T)
        inner.fields.B[1][I] = T(0.2) + (D >= 2 ? T(0.03) * cos(y) : zero(T))
        inner.fields.B[2][I] = T(0.1) * sin(x) + (D >= 3 ? T(0.02) * cos(z) : zero(T))
        inner.fields.B[3][I] = T(0.07) * cos(x) + (D >= 2 ? T(0.04) * sin(y) : zero(T))
    end
    st = CollisionalHybridStepper(
        inner,
        species;
        dt = T(0.01),
        rng = R(419),
        coulomb_log = T(strong ? 1e8 : 0.05),
        epsilon0 = one(T),
        weight_to_number = one(T),
    )
    init!(st, species)
    return st, species
end
