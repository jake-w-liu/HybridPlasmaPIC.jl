using HybridPlasmaPIC, Random, Test

# Both collisional wrappers must construct, initialize and advance through their
# public constructors; a type-header error makes construction impossible.
function collisional_wrapper_species(::Type{T}) where {T}
    g = FourierGrid((8,), (one(T),))
    electrons = ParticleSet{1,T}(4; q = -one(T), m = one(T))
    ions = ParticleSet{1,T}(4; q = one(T), m = T(4))
    for ps in (electrons, ions)
        ps.x[1] .= T[0.1, 0.35, 0.6, 0.85]
        ps.weight .= one(T)
    end
    electrons.v[1] .= T[0.01, -0.02, 0.015, -0.005]
    return g, [electrons, ions]
end

@testset "Collisional wrappers construct and advance" begin
    for T in (Float32, Float64), rng in (Xoshiro(1), MersenneTwister(1))
        g, species = collisional_wrapper_species(T)
        inner = EMPIC(g, species; mobile = Bool[true, true])
        em = CollisionalEMPICStepper(
            inner,
            species;
            dt = T(0.01),
            rng,
            coulomb_log = T(10),
            epsilon0 = one(T),
            weight_to_number = one(T),
        )
        @test em isa CollisionalEMPICStepper
        init!(em, species)
        step!(em, species, T(0.01))
        @test all(ps -> all(isfinite, ps.v[1]), species)

        g, species = collisional_wrapper_species(T)
        es = CollisionalElectrostaticPIC(
            ElectrostaticPIC(g, species),
            species;
            dt = T(0.01),
            rng,
            coulomb_log = T(10),
            epsilon0 = one(T),
            weight_to_number = one(T),
        )
        @test es isa CollisionalElectrostaticPIC
        init!(es, species)
        step!(es, species, T(0.01))
        @test all(ps -> all(isfinite, ps.v[1]), species)
    end
end
