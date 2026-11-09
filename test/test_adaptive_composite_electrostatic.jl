using Test
using LinearAlgebra
using HybridPlasmaPIC

# Small independent dense references are test-only. Production uses the tree and
# matrix-free operator; these references never call its assembly/cloud kernels.
function _ace_ref_constraints(mesh)
    D, R = length(mesh.L), length(mesh.vertices)
    vertices = [Float64.(v) for v in mesh.vertices]
    lookup = Dict(v => i for (i, v) in enumerate(vertices))
    bits = vec([Tuple(I) .- 1 for I in CartesianIndices(ntuple(_ -> 2, D))])
    C = zeros(R, R)
    hanging = Int[]
    for (v, x) in enumerate(vertices)
        incident = [
            k for k in eachindex(mesh.addresses) if
            all(d -> Float64(mesh.lower[k][d]) <= x[d] <= Float64(mesh.upper[k][d]), 1:D) &&
            any(d -> x[d] != Float64(mesh.lower[k][d]) && x[d] != Float64(mesh.upper[k][d]), 1:D)
        ]
        isempty(incident) && continue
        push!(hanging, v)
        k = first(sort(incident; by = j -> mesh.addresses[j]))
        lo, hi = Float64.(mesh.lower[k]), Float64.(mesh.upper[k])
        for bit in bits
            coordinate = ntuple(d -> bit[d] == 0 ? lo[d] : hi[d], D)
            C[v, lookup[coordinate]] += prod(
                d ->
                    bit[d] == 0 ? (hi[d] - x[d]) / (hi[d] - lo[d]) :
                    (x[d] - lo[d]) / (hi[d] - lo[d]),
                1:D,
            )
        end
    end
    masters = setdiff(collect(1:R), hanging)
    P = zeros(R, length(masters))
    P[masters, :] = Matrix{Float64}(I, length(masters), length(masters))
    isempty(hanging) || (
        P[hanging, :] =
            (Matrix{Float64}(I, length(hanging), length(hanging)) - C[hanging, hanging]) \
            C[hanging, masters]
    )
    free = [
        i for (i, v) in enumerate(masters) if all(d -> 0 < vertices[v][d] < Float64(mesh.L[d]), 1:D)
    ]
    return P, masters, free, bits, lookup
end
function _ace_ref_integral(x, a, lo, hi, f; knots = ())
    first, last = max(lo, x - a), min(hi, x + a)
    first < last || return 0.0
    breaks = sort(
        unique(vcat([first, last], [v for v in (x - a, x, x + a, knots...) if first < v < last])),
    )
    value = 0.0
    for i = 1:length(breaks)-1
        l, u = breaks[i], breaks[i+1]
        for sign in (-1, 1)
            y = (l + u) / 2 + sign * (u - l) / (2sqrt(3.0))
            density = max(0.0, 1 - abs((y - x) / a)) / a
            value += (u - l) / 2 * density * f(y)
        end
    end
    return value
end
function _ace_ref_moment(x, a, lo, hi, bit, L; tails = true)
    lower = tails && lo == 0 ? min(lo, x - a) : lo
    upper = tails && hi == L ? max(hi, x + a) : hi
    return _ace_ref_integral(
        x,
        a,
        lower,
        upper,
        y -> begin
            point = clamp(y, 0.0, L)
            bit == 0 ? (hi - point) / (hi - lo) : (point - lo) / (hi - lo)
        end;
        knots = (0.0, L, lo, hi),
    )
end
function _ace_ref_field(mesh, species, clouds, background)
    D = length(mesh.L)
    P, masters, free, bits, lookup = _ace_ref_constraints(mesh)
    K = zeros(length(masters), length(masters))
    rawsource = zeros(length(mesh.vertices))
    for k in eachindex(mesh.addresses)
        lo, hi = Float64.(mesh.lower[k]), Float64.(mesh.upper[k])
        widths = hi .- lo
        ids = [lookup[ntuple(d -> bit[d] == 0 ? lo[d] : hi[d], D)] for bit in bits]
        for quadrature in Iterators.product(ntuple(_ -> (-1, 1), D)...)
            xi = ntuple(d -> (1 + quadrature[d] / sqrt(3.0)) / 2, D)
            gradients = zeros(length(bits), D)
            for (a, bit) in enumerate(bits)
                N = prod(d -> bit[d] == 0 ? 1 - xi[d] : xi[d], 1:D)
                rawsource[ids[a]] += Float64(background) * prod(widths) * N / (2^D)
                for d = 1:D
                    gradients[a, d] =
                        (bit[d] == 0 ? -1 : 1) / widths[d] *
                        prod(j -> j == d ? 1.0 : (bit[j] == 0 ? 1 - xi[j] : xi[j]), 1:D)
                end
            end
            globalgrad = transpose(P[ids, :]) * gradients
            K .+= (prod(widths) / (2^D)) .* (globalgrad * transpose(globalgrad))
        end
        for ps in species, p in eachindex(ps.weight), (a, bit) in enumerate(bits)
            moment = prod(
                d -> _ace_ref_moment(
                    Float64(ps.x[d][p]),
                    Float64(clouds[d].half_width),
                    lo[d],
                    hi[d],
                    bit[d],
                    Float64(mesh.L[d]),
                ),
                1:D,
            )
            rawsource[ids[a]] += Float64(ps.q) * Float64(ps.weight[p]) * moment
        end
    end
    source = transpose(P) * rawsource
    potential = K[free, free] \ source[free]
    return (
        K = K,
        source = source,
        phi = potential,
        P = P,
        masters = masters,
        free = free,
        bits = bits,
        lookup = lookup,
    )
end
function _ace_ref_force(mesh, reference, phi, ps, p, clouds)
    D = length(mesh.L)
    u = zeros(length(reference.masters))
    u[reference.free] = Float64.(phi)
    raw = reference.P * u
    E = zeros(3)
    for k in eachindex(mesh.addresses)
        lo, hi = Float64.(mesh.lower[k]), Float64.(mesh.upper[k])
        for bit in reference.bits
            id = reference.lookup[ntuple(d -> bit[d] == 0 ? lo[d] : hi[d], D)]
            for d = 1:D
                mass = _ace_ref_integral(
                    Float64(ps.x[d][p]),
                    Float64(clouds[d].half_width),
                    lo[d],
                    hi[d],
                    _ -> 1.0,
                )
                derivative = (bit[d] == 0 ? -1 : 1) * mass / (hi[d] - lo[d])
                for j = 1:D
                    j == d && continue
                    derivative *= _ace_ref_moment(
                        Float64(ps.x[j][p]),
                        Float64(clouds[j].half_width),
                        lo[j],
                        hi[j],
                        bit[j],
                        Float64(mesh.L[j]);
                        tails = false,
                    )
                end
                E[d] -= raw[id] * derivative
            end
        end
    end
    return E
end
_ace_test_value(x) = x
_ace_test_value(x::Tuple) = map(_ace_test_value, x)
_ace_test_value(x::NamedTuple) = NamedTuple{keys(x)}(map(_ace_test_value, values(x)))
_ace_test_value(x::AbstractArray) = map(_ace_test_value, x)
function _ace_test_value(
    x::Union{
        AdaptiveCompositeESSolveReport,
        AdaptiveCompositeESStepReport,
        AdaptiveCompositeESRegridReport,
        HybridPlasmaPIC._AdaptiveCompositeCrossing,
    },
)
    return NamedTuple{fieldnames(typeof(x))}(
        ntuple(i -> _ace_test_value(getfield(x, i)), fieldcount(typeof(x))),
    )
end
_ace_test_snapshot(st, species) = _ace_test_value(HybridPlasmaPIC._checkpoint_state(st, species))
function _ace_test_population(::Val{D}, ::Type{T}; first_count = 1) where {D,T}
    a, b, tracer, empty = ParticleSet{D,T}(first_count; q = T(1), m = T(2)),
    ParticleSet{D,T}(1; q = T(-0.75), m = T(1.5)),
    ParticleSet{D,T}(1; q = T(0.5), m = T(1)),
    ParticleSet{D,T}(0)
    for p in eachindex(a.weight), d = 1:D
        a.x[d][p] = T((2.75, 3.25, 2.5)[d])
        a.v[d][p] = T((0.125, 0.1875, 0.0625)[d])
        a.weight[p] = T(1.25)
    end
    for d = 1:D
        b.x[d][1], b.v[d][1] = T((0.25, 0.375, 0.5)[d]), T(0.03125)
        tracer.x[d][1], tracer.v[d][1] = T((1.5, 4.75, 2.25)[d]), T(0.0625)
    end
    b.weight[1], tracer.weight[1] = T(0.625), zero(T)
    D == 2 && fill!(a.v[3], T(0.25))
    return [a, b, tracer, empty]
end
@testset "Adaptive Cartesian composite electrostatic essentials" begin
    H = HybridPlasmaPIC
    for D in (2, 3), T in (Float32, Float64)
        @testset "D=$D $T" begin
            base = CartesianAMRHierarchy(ntuple(_ -> 3, D), ntuple(_ -> T(6), D))
            fine = CartesianAMRHierarchy(base.nbase, base.L, ((0, ntuple(_ -> 2, D)),))
            @test length(fine.addresses) == 3^D + 2^D - 1
            @test sum(fine.volumes) == 6^D
            @test_throws ArgumentError CartesianAMRHierarchy(
                base.nbase,
                base.L,
                ((0, ntuple(_ -> 2, D)), (1, ntuple(_ -> 3, D))),
            )
            if D == 3 && T === Float64
                specs = Tuple{Int,NTuple{D,Int}}[]
                for I in CartesianIndices(base.nbase)
                    push!(specs, (0, Tuple(I)))
                end
                push!(specs, (1, ntuple(_ -> 3, D)))
                deeper = CartesianAMRHierarchy(base.nbase, base.L, specs)
                @test maximum(first, deeper.addresses) == 2
                @test sum(deeper.volumes) == 6^D
                @test all(row -> sum(last, row) == 1, deeper.rows)
            end
            clouds = ntuple(
                d ->
                    d == 2 ? YeeTriangularCloud(Float32(0.5)) :
                    YeeTriangularCloud(Float64(d == 1 ? 0.625 : 0.375)),
                D,
            )
            species = _ace_test_population(Val(D), T)
            st = AdaptiveCompositeElectrostaticPIC(
                fine,
                species;
                clouds,
                background_charge = T(0.03125),
                rtol = T(16) * eps(T),
            )
            init!(st, species)
            reference = _ace_ref_field(fine, species, clouds, st.background_charge)
            @test fine.masters == reference.masters && fine.free == reference.free
            actualP = zeros(length(fine.vertices), length(fine.masters))
            for i in eachindex(fine.rows), (j, c) in fine.rows[i]
                actualP[i, j] = Float64(c)
            end
            @test isapprox(actualP, reference.P; rtol = 1e-13, atol = 1e-13)
            @test isapprox(reference.K, transpose(reference.K); rtol = 1e-13, atol = 1e-13)
            @test maximum(abs, reference.K * ones(length(fine.masters))) < 1e-12
            @test minimum(eigvals(Symmetric(reference.K[reference.free, reference.free]))) > 0
            diagnostic = composite_electrostatic_diagnostics(st, species)
            @test isapprox(
                Float64.(diagnostic.nodal_source),
                reference.source;
                rtol = 2e-12,
                atol = 2e-12,
            )
            @test isapprox(Float64.(st.phi), reference.phi; rtol = 4096eps(T), atol = 4096eps(T))
            @test diagnostic.inventory.particle == 1 * T(1.25) - T(0.75) * T(0.625)
            @test !iszero(diagnostic.inventory.tail)
            @test diagnostic.inventory.inside + diagnostic.inventory.tail ==
                  diagnostic.inventory.particle
            @test iszero(diagnostic.accounting_defect)
            u = zeros(length(fine.masters))
            u[fine.free] = Float64.(st.phi)
            @test isapprox(
                Float64(diagnostic.field_energy),
                dot(u, reference.K * u) / 2;
                rtol = 2e-12,
                atol = 2e-12,
            )
            raw = H._acq_raw_potential(fine, st.phi)
            forces = [
                _ace_ref_force(fine, reference, st.phi, ps, p, clouds) for ps in species for
                p in eachindex(ps.weight)
            ]
            for (row, ps) in enumerate(species), p in eachindex(ps.weight)
                expected = _ace_ref_force(fine, reference, st.phi, ps, p, clouds)
                actual =
                    H._acc_gather(fine, raw, ntuple(d -> H._exact_rational(ps.x[d][p]), D), clouds)
                @test isapprox(Float64.(collect(actual)), expected; rtol = 2e-11, atol = 2e-11)
            end
            @test all(isfinite, composite_electric_field(st, ntuple(_ -> zero(T), D)))
            if D == 2 && T === Float64
                origin = st.initial_energy
                regrid!(st, species, base)
                regrid!(st, species, fine)
                @test st.step == 0 && st._history === nothing && st.initial_energy == origin
                @test iszero(st.physical_defect) &&
                      iszero(composite_electrostatic_diagnostics(st, species).accounting_defect)
                wideclouds = ntuple(_ -> YeeTriangularCloud(Float64(9)), D)
                wideparticles = _ace_test_population(Val(D), T)
                wide = AdaptiveCompositeElectrostaticPIC(
                    base,
                    wideparticles;
                    clouds = wideclouds,
                    rtol = T(16) * eps(T),
                )
                init!(wide, wideparticles)
                wideref = _ace_ref_field(base, wideparticles, wideclouds, wide.background_charge)
                @test isapprox(
                    Float64.(composite_electrostatic_diagnostics(wide, wideparticles).nodal_source),
                    wideref.source;
                    rtol = 2e-12,
                    atol = 2e-12,
                )
            end
            before = _ace_test_snapshot(st, species)
            mktempdir() do dir
                path = joinpath(dir, "initialized.run")
                archive_run(path, st, species; rng_seed = 0)
                load_archive!(st, species, path)
                @test isequal(_ace_test_snapshot(st, species), before)
            end
            @test step!(st, species, -zero(T)) === st
            @test isequal(_ace_test_snapshot(st, species), before)
            @test regrid!(st, species, fine) === st
            @test isequal(_ace_test_snapshot(st, species), before)
            h = T(1 / 32)
            oldx, oldv = [map(copy, ps.x) for ps in species], [map(copy, ps.v) for ps in species]
            vhalf = NTuple{3,T}[]
            cursor = 0
            for ps in species, p in eachindex(ps.weight)
                cursor += 1
                push!(
                    vhalf,
                    ntuple(
                        c -> T(
                            Float64(ps.v[c][p]) +
                            Float64(h) / 2 * Float64(ps.q) / Float64(ps.m) * forces[cursor][c],
                        ),
                        3,
                    ),
                )
            end
            step!(st, species, h)
            endref = _ace_ref_field(fine, species, clouds, st.background_charge)
            @test isapprox(Float64.(st.phi), endref.phi; rtol = 4096eps(T), atol = 4096eps(T))
            cursor = 0
            for (s, ps) in enumerate(species), p in eachindex(ps.weight)
                cursor += 1
                endforce = _ace_ref_force(fine, endref, st.phi, ps, p, clouds)
                for d = 1:D
                    @test isapprox(
                        ps.x[d][p],
                        T(oldx[s][d][p] + h * vhalf[cursor][d]);
                        rtol = 64eps(T),
                        atol = 64eps(T),
                    )
                    @test ps.x[d][p] != oldx[s][d][p]
                end
                for c = 1:3
                    expected = T(
                        Float64(vhalf[cursor][c]) +
                        Float64(h) / 2 * Float64(ps.q) / Float64(ps.m) * endforce[c],
                    )
                    @test isapprox(ps.v[c][p], expected; rtol = 64eps(T), atol = 64eps(T))
                end
            end
            if D == 2
                @test species[1].v[3][1] == oldv[1][3][1]
            end
            r = st._state.step_report
            @test r.total_defect == r.field_defect + r.kick_defect + r.trajectory_defect
            @test !iszero(r.trajectory_defect)
            @test_throws ArgumentError step!(st, species, T(1000))
            stepped = _ace_test_snapshot(st, species)
            mktempdir() do dir
                path = joinpath(dir, "physical-step.chk")
                save_checkpoint(path, st, species)
                load_checkpoint!(st, species, path)
                @test isequal(_ace_test_snapshot(st, species), stepped)
            end
            @test_throws ArgumentError step!(st, species, -h)
            @test isequal(_ace_test_snapshot(st, species), stepped)
            regrid!(st, species, base)
            @test !iszero(st.last_regrid_report.mesh_work)
            history = _ace_test_value(st._history)
            regrid!(st, species, fine)
            @test isequal(_ace_test_value(st._history), history)
            @test st._history.state.refinements == fine.refinements
            @test iszero(composite_electrostatic_diagnostics(st, species).accounting_defect)
            saved = _ace_test_snapshot(st, species)
            original = H._checkpoint_state(st, species)
            aliased = merge(original, (current = merge(original.current, (phi = st.phi,)),))
            @test_throws ArgumentError H._restore_checkpoint_state!(st, species, aliased)
            @test_throws ArgumentError H._restore_checkpoint_state!(
                st,
                species,
                merge(original, (mesh_work = original.mesh_work + 1,)),
            )
            @test isequal(_ace_test_snapshot(st, species), saved)
            mktempdir() do dir
                rawpath, archivepath = joinpath(dir, "amr.chk"), joinpath(dir, "amr.run")
                save_checkpoint(rawpath, st, species)
                archive_run(archivepath, st, species; rng_seed = 0)
                other = _ace_test_population(Val(D), T; first_count = 2)
                destination = AdaptiveCompositeElectrostaticPIC(
                    base,
                    other;
                    clouds,
                    background_charge = T(0.03125),
                    rtol = T(16) * eps(T),
                )
                init!(destination, other)
                for ps in other
                    foreach(a -> fill!(a, T(NaN)), ps.x)
                    fill!(ps.weight, T(NaN))
                    ps.q, ps.m = T(NaN), -one(T)
                end
                fill!(destination.phi, T(NaN))
                fill!(destination._state.phi, T(NaN))
                destination.time, destination.step, destination.last_dt = T(Inf), -1, T(Inf)
                owner = other[1]
                load_checkpoint!(destination, other, rawpath)
                @test other[1] === owner && length(other[1].weight) == 1
                @test isequal(_ace_test_snapshot(destination, other), saved)
                step!(st, species, h)
                step!(destination, other, h)
                @test isequal(
                    _ace_test_snapshot(destination, other),
                    _ace_test_snapshot(st, species),
                )
                @test destination._history === nothing && destination.step == 2
                load_archive!(destination, other, archivepath)
                @test isequal(_ace_test_snapshot(destination, other), saved)
                if D == 2 && T === Float64
                    adapt!(
                        destination,
                        other;
                        refine_threshold = zero(T),
                        coarsen_threshold = zero(T),
                        max_level = 2,
                        buffer = 0,
                    )
                    @test maximum(first, destination.mesh.addresses) == 2
                end
                adapt!(
                    destination,
                    other;
                    refine_threshold = floatmax(T),
                    coarsen_threshold = floatmax(T) / 2,
                    max_level = 2,
                    buffer = 0,
                )
                @test isempty(destination.mesh.refinements)
                @test iszero(
                    composite_electrostatic_diagnostics(destination, other).accounting_defect,
                )
            end
            # Simultaneous internal face/edge/corner crossing has one time event.
            neutral = ParticleSet{D,T}(1; q = zero(T))
            for d = 1:D
                neutral.x[d][1], neutral.v[d][1] = T(2.75), one(T)
            end
            crossing =
                AdaptiveCompositeElectrostaticPIC(fine, [neutral]; clouds, rtol = T(16) * eps(T))
            init!(crossing, [neutral])
            step!(crossing, [neutral], T(0.5))
            @test length(crossing._state.step_report.crossings) == 1
            event = only(crossing._state.step_report.crossings)
            @test event.parameter == 1 // 2 && all(event.axes)
            for d = 1:D
                neutral.x[d][1], neutral.v[d][1] = T(2.5), one(T)
            end
            init!(crossing, [neutral])
            step!(crossing, [neutral], T(0.5))
            @test only(crossing._state.step_report.crossings).parameter == 1
            for d = 1:D
                neutral.x[d][1], neutral.v[d][1] = T(3), -one(T)
            end
            init!(crossing, [neutral])
            step!(crossing, [neutral], T(0.5))
            @test only(crossing._state.step_report.crossings).parameter == 0
            # A true staged physical-wall failure leaves a usable initial state.
            fast = _ace_test_population(Val(D), T)
            fast[1].v[1][1] = T(2000)
            wall = AdaptiveCompositeElectrostaticPIC(fine, fast; clouds, rtol = T(16) * eps(T))
            init!(wall, fast)
            wallbefore = _ace_test_snapshot(wall, fast)
            @test_throws ArgumentError step!(wall, fast, T(1 / 32))
            @test isequal(_ace_test_snapshot(wall, fast), wallbefore)
            step!(wall, fast, T(1 / 65536))
            @test wall.step == 1
            # A real failed endpoint/regrid solve: one coarse free unknown solves
            # in1 CG iteration, while the refined nonuniform source needs more.
            limitedmesh = CartesianAMRHierarchy(ntuple(_ -> 2, D), ntuple(_ -> T(4), D))
            q = ParticleSet{D,T}(1)
            for d = 1:D
                q.x[d][1] = T((0.75, 1.125, 0.625)[d])
            end
            limited = AdaptiveCompositeElectrostaticPIC(
                limitedmesh,
                [q];
                clouds,
                rtol = T(16) * eps(T),
                maxiter = 1,
            )
            init!(limited, [q])
            limitedbefore = _ace_test_snapshot(limited, [q])
            proposed =
                CartesianAMRHierarchy(limitedmesh.nbase, limitedmesh.L, ((0, ntuple(_ -> 1, D)),))
            @test_throws ArgumentError regrid!(limited, [q], proposed)
            @test isequal(_ace_test_snapshot(limited, [q]), limitedbefore)
            step!(limited, [q], T(1 / 128))
            @test limited.step == 1
        end
    end
end
