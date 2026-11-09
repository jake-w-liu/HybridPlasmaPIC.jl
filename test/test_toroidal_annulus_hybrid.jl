module ToroidalAnnulusHybridTests
using Test
using HybridPlasmaPIC
const H = HybridPlasmaPIC
const Q = Rational{BigInt}

function fixture(::Type{T}) where {T}
    mesh = ToroidalAnnulusGrid((3, 3, 3); R0 = 2, r_inner = 1 / 4, r_outer = 1 / 2, T = T)
    count = (mesh.n[1] - 1) * mesh.n[2] * mesh.n[3]
    species = [
        ParticleSet{3,T}(count; q = 1, m = 1),
        ParticleSet{3,T}(count; q = 2, m = 3),
        ParticleSet{3,T}(0; q = 3, m = 5),
    ]
    for (index, p) in enumerate(species[1:2])
        load_toroidal_annulus_cold!(
            p,
            mesh;
            number_density = T(index == 1 ? 3 / 4 : 1 / 4),
            velocity = (T(1 / 64), -T(1 / 128), T(1 / 32)),
        )
        p.id .= UInt64.(index * count .+ (1:count))
        p.tag .= UInt32(index)
    end
    beta = toroidal_annulus_vacuum_flux(mesh; strength = T(1 / 2))
    pressure = zeros(T, mesh.n)
    for I in CartesianIndices(pressure)
        u, v = mesh.op.coordinates[2][I[2]], mesh.op.coordinates[3][I[3]]
        # Sole beta3 stays phi-independent so its initial logical divergence is zero.
        beta[3][I] *= one(T) + cospi(T(2) * u) / T(16)
        pressure[I] = T(1 / 8) * (one(T) + sinpi(T(2) * u) / T(8) + cospi(T(2) * v) / T(8))
    end
    st = ToroidalAnnulusHybridPIC(
        mesh,
        species;
        beta0 = beta,
        pe0 = pressure,
        pressure_reservoir = (T(1 / 4), T(1 / 32)),
    )
    return mesh, st, species
end

# Independent eight-corner CIC sum, with the declared stored coordinate/weight
# contract. Only final IEEE rational rounding uses the shared scalar converter.
function density_reference(mesh, species, ::Type{T}) where {T}
    charge = fill(zero(Q), mesh.n)
    for p in species, marker in eachindex(p.weight)
        s = ntuple(d -> (d == 1 ? p.x[d][marker] : mod(p.x[d][marker], one(T))) / mesh.dx[d], 3)
        base = map(x -> floor(Int, x), s)
        fraction = ntuple(d -> s[d] - base[d], 3)
        for corner in CartesianIndices((2, 2, 2))
            weight = Q(p.q) * Q(p.weight[marker])
            for d = 1:3
                weight *= Q(corner[d] == 1 ? one(T) - fraction[d] : fraction[d])
            end
            index = ntuple(
                d ->
                    d == 1 ? clamp(base[d] + corner[d], 1, mesh.n[d]) :
                    mod(base[d] + corner[d] - 1, mesh.n[d]) + 1,
                3,
            )
            charge[index...] += weight
        end
    end
    density = similar(charge, T)
    for I in CartesianIndices(charge)
        volume =
            Q(mesh.tau)^2 *
            Q(mesh.r[I[1]]) *
            Q(mesh.R[I[1], I[2]]) *
            Q(mesh.op.sbp.H[I[1]]) *
            Q(mesh.dx[2]) *
            Q(mesh.dx[3])
        density[I] = H._from_exact(T, charge[I] / volume)
    end
    return density
end

function leaves(st, species)
    particles = Tuple(a for p in species for a in (p.x..., p.v..., p.weight, p.id, p.tag))
    return (st.beta..., H._tah_field_arrays(st.fields)..., particles...)
end

function make_stale!(st, species)
    foreach(a -> fill!(a, NaN), (st.beta..., H._tah_field_arrays(st.fields)...))
    for p in species
        foreach(a -> fill!(a, NaN), (p.x..., p.v..., p.weight))
        p.q, p.m = NaN, NaN
    end
    st.initialized, st.time, st.step, st.last_dt = false, NaN, -3, NaN
    return nothing
end

@testset "Coupled toroidal annulus essentials" begin
    for T in (Float32, Float64)
        @testset "$T" begin
            mesh, st, species = fixture(T)
            point = (T(1 / 8), zero(T), zero(T))
            @test to_cartesian(mesh, point) == (T(19 / 8), zero(T), zero(T))
            @test to_logical(mesh, to_cartesian(mesh, point)) == point
            @test_throws ArgumentError ToroidalAnnulusGrid(
                (3, 4, 3);
                R0 = 2,
                r_inner = 1 / 4,
                r_outer = 1 / 2,
                T = T,
            )
            seed = toroidal_annulus_vacuum_flux(mesh; strength = T(1 / 2))
            @test all(iszero, seed[1]) && all(iszero, seed[2])
            @test all(
                I ->
                    seed[3][I] ==
                    H._from_exact(T, Q(mesh.tau) * Q(mesh.r[I[1]]) / (2Q(mesh.R[I[1], I[2]]))),
                CartesianIndices(mesh.n),
            )
            # The first quiet marker is at the midpoint of the first actual cell.
            p = first(species)
            marker = ntuple(
                d -> H._from_exact(
                    T,
                    (Q(mesh.op.coordinates[d][1]) + Q(mesh.op.coordinates[d][2])) / 2,
                ),
                3,
            )
            @test ntuple(d -> p.x[d][1], 3) == marker
            r = H._from_exact(T, Q(mesh.r_inner) + Q(marker[1]))
            major = H._from_exact(T, Q(mesh.R0) + Q(r) * Q(cospi(T(2) * marker[2])))
            cell = prod(Q(mesh.op.coordinates[d][2]) - Q(mesh.op.coordinates[d][1]) for d = 1:3)
            @test p.weight[1] ==
                  H._from_exact(T, Q(T(3 / 4)) * Q(mesh.tau)^2 * Q(r) * Q(major) * cell)
            invalid = ParticleSet{3,T}(length(p.weight); q = 0, m = 1)
            foreach(a -> fill!(a, NaN), (invalid.x..., invalid.v..., invalid.weight))
            unchanged =
                map(copy, (invalid.x..., invalid.v..., invalid.weight, invalid.id, invalid.tag))
            @test_throws ArgumentError load_toroidal_annulus_cold!(
                invalid,
                mesh;
                number_density = 1,
            )
            @test isequal(
                (invalid.x..., invalid.v..., invalid.weight, invalid.id, invalid.tag),
                unchanged,
            )
            expected_density = density_reference(mesh, species, T)
            identities = leaves(st, species)
            @test init!(st, species) === st
            @test isequal(st.fields.n, expected_density)
            @test all(>(zero(T)), st.fields.n)
            @test st.last_endpoint.thermal_controller[1] > 0
            @test st.last_endpoint.thermal_controller[2] < 0
            initial = H._checkpoint_state(st, species)
            # A syntactically valid zero-step receipt must still describe the
            # actual contact, even before any completed interval exists.
            endpoint = st.last_endpoint
            st.last_endpoint = typeof(endpoint)(
                ntuple(
                    i ->
                        fieldnames(typeof(endpoint))[i] === :nu ? (zero(Q), zero(Q)) :
                        getfield(endpoint, i),
                    fieldcount(typeof(endpoint)),
                )...,
            )
            @test_throws ArgumentError step!(st, species, zero(T))
            @test isequal(st.fields, initial.fields)
            st.last_endpoint = endpoint
            @test isequal(H._checkpoint_state(st, species), initial)
            @test step!(st, species, -zero(T)) === st
            @test isequal(H._checkpoint_state(st, species), initial)
            @test_throws ArgumentError step!(st, species, -one(T))
            @test isequal(H._checkpoint_state(st, species), initial)
            # Deliberately exceed a measured incoming SAT frequency, reaching
            # the private coupled RHS rather than changing the public fixture.
            frequency = maximum(st.last_endpoint.nu)
            @test frequency > 0
            rejected = try
                step!(st, species, T(4 / frequency))
                nothing
            catch err
                err
            end
            @test rejected isa ArgumentError && occursin("SAT-only", sprint(showerror, rejected))
            @test isequal(H._checkpoint_state(st, species), initial)
            h = T(1 / 1024)
            @test step!(st, species, h) === st
            @test st.time == h && st.step == 1 && st.last_dt == h
            @test any(i -> !isequal(species[i].v, initial.species[i].v), 1:2)
            @test !isequal(st.beta, initial.beta)
            @test !isequal(st.fields.pe, initial.fields.pe)
            @test all(pair -> pair[1] === pair[2], zip(identities, leaves(st, species)))
            @test all(
                side -> all(
                    pair -> Q(pair[1]) == Q(pair[2]),
                    zip(st.beta[1][side, :, :], initial.beta[1][side, :, :]),
                ),
                (1, mesh.n[1]),
            )
            d = toroidal_annulus_diagnostics(st)
            @test iszero(d.accounting_defect)
            @test d.energy.total - d.initial_energy + d.cumulative_outward_work -
                  d.cumulative_controller_work == d.cumulative_stage_defect + d.cumulative_rk_defect
            @test all(r -> r.divergence <= r.divergence_tolerance, st.last_report.stages)
            @test step!(st, species, 3h / 2) === st
            saved = H._checkpoint_state(st, species)
            mktempdir() do directory
                raw, archive =
                    joinpath(directory, "annulus.bin"), joinpath(directory, "annulus.run")
                save_checkpoint(raw, st, species)
                archive_run(archive, st, species; rng_seed = 0)
                step!(st, species, 3h / 4)
                future = H._checkpoint_state(st, species)
                for format in (:raw, :archive)
                    _, restored, populations = fixture(T)
                    init!(restored, populations)
                    retained = leaves(restored, populations)
                    make_stale!(restored, populations)
                    result =
                        format === :raw ? load_checkpoint!(restored, populations, raw) :
                        load_archive!(restored, populations, archive)
                    @test result === restored
                    @test isequal(H._checkpoint_state(restored, populations), saved)
                    @test all(
                        pair -> pair[1] === pair[2],
                        zip(retained, leaves(restored, populations)),
                    )
                    step!(restored, populations, 3h / 4)
                    @test isequal(H._checkpoint_state(restored, populations), future)
                end
                @test_throws ArgumentError H._restore_checkpoint_state!(
                    st,
                    species,
                    merge(future, (schema = 2,)),
                )
                @test isequal(H._checkpoint_state(st, species), future)
                @test_throws ArgumentError H._restore_checkpoint_state!(
                    st,
                    species,
                    merge(future, (beta = st.beta,)),
                )
                @test isequal(H._checkpoint_state(st, species), future)
            end
        end
    end
end
end
