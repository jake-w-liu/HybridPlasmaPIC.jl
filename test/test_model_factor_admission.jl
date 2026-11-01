module ModelFactorAdmissionTests
using HybridPlasmaPIC, Test
const H = HybridPlasmaPIC
const Q = Rational{BigInt}

function factors(::Type{T}, overflow) where {T}
    tiny = nextfloat(zero(T))
    if overflow
        return (q = one(T), m = tiny, dt = 2 * floatmin(T), E = ldexp(one(T), -precision(T) - 8))
    end
    k, b = T === Float32 ? (80, 60) : (600, 400)
    return (q = tiny, m = ldexp(one(T), k), dt = ldexp(one(T), b), E = ldexp(one(T), b))
end

function particles(::Type{T}, q, m) where {T}
    ps = ParticleSet{1,T}(4; q, m)
    ps.x[1] .= T[1, 3, 5, 7] ./ T(8)
    return ps
end

function check_first_half_step(ps, x0, expected)
    @test ps.x == x0 # The imposed force is transverse to this 1D grid.
    @test all(iszero, ps.v[1]) && all(iszero, ps.v[2])
    @test all(==(expected), ps.v[3])
end

@testset "Hybrid model admission retains original charge and mass" begin
    for T in (Float32, Float64), St in (HybridStepper, CAMCLStepper), overflow in (true, false)
        f = factors(T, overflow)
        @test overflow ? isinf(f.q / f.m) : iszero(f.q / f.m)
        # The scalar interface keeps its proton-normalized moments; the
        # collection interface supports both charge magnitudes.
        for collection in (overflow ? (false, true) : (true,))
            g = FourierGrid((4,), (one(T),))
            ps = particles(T, f.q, f.m)
            populations = collection ? [ps] : ps
            st = St(g, HybridModel(IsothermalElectrons(0)), CIC(), 0)
            initialize = St === HybridStepper ? init! : init_camcl!
            advance = St === HybridStepper ? step! : step_camcl!
            @test initialize(st, populations) === st
            st.fields.E[3] .= f.E
            x0 = deepcopy(ps.x)
            # Constant E and B=0 give v^(1/2) = (dt/2)*(q/m)*E.
            exact = Q(f.dt) * Q(f.q) * Q(f.E) / (2 * Q(f.m))
            expected = T(exact)
            @test !iszero(expected) && Q(expected) == exact
            advance(st, populations, f.dt; NB = 2)
            check_first_half_step(ps, x0, expected)
            @test st.time[] == f.dt && st.step[] == 1 && st.last_dt[] == f.dt
            saved = deepcopy(H._checkpoint_state(st, populations))
            advance(st, populations, 2 * f.dt; NB = 2)
            continued = deepcopy(H._checkpoint_state(st, populations))
            @test all(==(expected), ps.v[3]) # The recomputed self-field is zero.
            H._restore_checkpoint_state!(st, populations, saved)
            advance(st, populations, 2 * f.dt; NB = 2)
            @test H._checkpoint_state(st, populations) == continued
        end
    end
end

@testset "Electrostatic collection admission retains original factors" begin
    for T in (Float32, Float64), overflow in (true, false), sign in (-1, 1)
        f = factors(T, overflow)
        q = sign * f.q
        @test overflow ? isinf(q / f.m) : iszero(q / f.m)
        g = FourierGrid((4,), (one(T),))
        ps = particles(T, q, f.m)
        species = [ps]
        es = ElectrostaticPIC(g, species; neutrality = :neutralize)
        @test init_espic!(es, species) === es
        es.E[3] .= f.E
        x0 = deepcopy(ps.x)
        exact = Q(f.dt) * Q(q) * Q(f.E) / (2 * Q(f.m))
        expected = T(exact)
        @test !iszero(expected) && Q(expected) == exact
        step_espic!(es, species, f.dt)
        check_first_half_step(ps, x0, expected)
        @test es.primed[] && es.last_dt[] == f.dt
        @test electrostatic_diagnostics(es, species).particle_charge == 4 * Q(q)
        saved = deepcopy(H._checkpoint_state(es, species))
        step_espic!(es, species, 2 * f.dt)
        continued = deepcopy(H._checkpoint_state(es, species))
        @test all(==(expected), ps.v[3])
        H._restore_checkpoint_state!(es, species, saved)
        step_espic!(es, species, 2 * f.dt)
        @test H._checkpoint_state(es, species) == continued
        # Invalid original operands must still fail before publication.
        ps.m = zero(T)
        before = deepcopy((ps.x, ps.v, es.E, es.ne, es.primed[], es.last_dt[]))
        @test_throws ArgumentError step_espic!(es, species, f.dt)
        @test (ps.x, ps.v, es.E, es.ne, es.primed[], es.last_dt[]) == before
    end
end
end
