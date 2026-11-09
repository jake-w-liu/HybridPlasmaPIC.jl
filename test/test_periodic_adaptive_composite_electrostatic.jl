using Test
using LinearAlgebra
using HybridPlasmaPIC

# Tiny independent torus reference: quotient coordinates, dense trace equations,
# then physical gradient/density quadrature. No production P or kernel is used.
function _periodic_amr_ref_space(mesh)
    D = length(mesh.L)
    Q = Rational{BigInt}
    canonical(x) = ntuple(d -> mod(x[d],Q(mesh.L[d])),D)
    classes = sort!(unique(canonical.(mesh.vertices)))
    class_index = Dict(x => i for (i,x) in enumerate(classes))
    raw_class = [class_index[canonical(x)] for x in mesh.vertices]
    raw_index = Dict(x => i for (i,x) in enumerate(mesh.vertices))
    bits = vec([Tuple(I) .- 1 for I in CartesianIndices(ntuple(_ -> 2,D))])
    C = zeros(length(classes),length(classes))
    hanging = Int[]
    for (c,x) in enumerate(classes)
        candidates = Tuple{Int,NTuple{D,Q}}[]
        choices = ntuple(d -> iszero(x[d]) ? (zero(Q),Q(mesh.L[d])) : (x[d],),D)
        for image in Iterators.product(choices...),k in eachindex(mesh.addresses)
            all(d -> mesh.lower[k][d] <= image[d] <= mesh.upper[k][d],1:D) || continue
            any(d -> image[d] != mesh.lower[k][d] && image[d] != mesh.upper[k][d],1:D) && push!(candidates,(k,Tuple(image)))
        end
        isempty(candidates) && continue
        push!(hanging,c)
        k,image = first(sort!(candidates;by = pair -> (mesh.addresses[pair[1]],pair[2])))
        for bit in bits
            corner = ntuple(d -> bit[d] == 0 ? mesh.lower[k][d] : mesh.upper[k][d],D)
            coefficient = prod(d -> bit[d] == 0 ? (mesh.upper[k][d]-image[d])/(mesh.upper[k][d]-mesh.lower[k][d]) : (image[d]-mesh.lower[k][d])/(mesh.upper[k][d]-mesh.lower[k][d]),1:D)
            C[c,class_index[canonical(corner)]] += Float64(coefficient)
        end
    end
    master_classes = setdiff(collect(eachindex(classes)),hanging)
    Pc = zeros(length(classes),length(master_classes))
    Pc[master_classes,:] = Matrix{Float64}(I,length(master_classes),length(master_classes))
    isempty(hanging) || (Pc[hanging,:] = (Matrix{Float64}(I,length(hanging),length(hanging))-C[hanging,hanging])\C[hanging,master_classes])
    representatives = [findfirst(==(c),raw_class)::Int for c in master_classes]
    gauge = findfirst(c -> all(iszero,classes[c]),master_classes)::Int
    free = setdiff(collect(eachindex(master_classes)),[gauge])
    return (P = Pc[raw_class,:],masters = representatives,free = free,gauge = gauge,bits = bits,raw_index = raw_index)
end
function _periodic_amr_ref_moment(x,a,lo,hi,L,bit)
    value = 0.0
    images = ceil(Int,a/L)+2
    for image = -images:images
        centre = x+image*L
        left,right = max(lo,centre-a),min(hi,centre+a)
        left < right || continue
        knots = sort!(unique!([left,right,clamp(centre,left,right)]))
        for j = 1:length(knots)-1,sign in (-1,1)
            l,u = knots[j],knots[j+1]
            point = (l+u)/2+sign*(u-l)/(2sqrt(3.0))
            density = max(0.0,1-abs((point-centre)/a))/a
            basis = bit === nothing ? 1.0 : (bit == 0 ? (hi-point)/(hi-lo) : (point-lo)/(hi-lo))
            value += (u-l)/2*density*basis
        end
    end
    return value
end
function _periodic_amr_ref_field(mesh,species,clouds,density)
    D = length(mesh.L)
    space = _periodic_amr_ref_space(mesh)
    K = zeros(length(space.masters),length(space.masters))
    source = zeros(length(mesh.vertices))
    for k in eachindex(mesh.addresses)
        lo,hi = Float64.(mesh.lower[k]),Float64.(mesh.upper[k])
        h = hi .- lo
        ids = [space.raw_index[ntuple(d -> bit[d] == 0 ? mesh.lower[k][d] : mesh.upper[k][d],D)] for bit in space.bits]
        for point in Iterators.product(ntuple(_ -> (-1,1),D)...)
            xi = ntuple(d -> (1+point[d]/sqrt(3.0))/2,D)
            gradients = zeros(length(space.bits),D)
            for (a,bit) in enumerate(space.bits)
                source[ids[a]] += Float64(density)*prod(h)/(2^D)*prod(d -> bit[d] == 0 ? 1-xi[d] : xi[d],1:D)
                for d = 1:D
                    gradients[a,d] = (bit[d] == 0 ? -1 : 1)/h[d]*prod(j -> j == d ? 1.0 : (bit[j] == 0 ? 1-xi[j] : xi[j]),1:D)
                end
            end
            global_gradient = transpose(space.P[ids,:])*gradients
            K .+= prod(h)/(2^D).*(global_gradient*transpose(global_gradient))
        end
        for ps in species,p in eachindex(ps.weight),(a,bit) in enumerate(space.bits)
            moment = prod(d -> _periodic_amr_ref_moment(Float64(ps.x[d][p]),Float64(clouds[d].half_width),lo[d],hi[d],Float64(mesh.L[d]),bit[d]),1:D)
            source[ids[a]] += Float64(ps.q)*Float64(ps.weight[p])*moment
        end
    end
    b = transpose(space.P)*source
    phi = K[space.free,space.free]\b[space.free]
    return merge(space,(K = K,source = b,phi = phi))
end
function _periodic_amr_ref_force(mesh,reference,phi,ps,p,clouds)
    D = length(mesh.L)
    u = zeros(length(reference.masters));u[reference.free] = Float64.(phi)
    raw = reference.P*u
    field = zeros(3)
    for k in eachindex(mesh.addresses),bit in reference.bits
        lo,hi = Float64.(mesh.lower[k]),Float64.(mesh.upper[k])
        raw_index = reference.raw_index[ntuple(d -> bit[d] == 0 ? mesh.lower[k][d] : mesh.upper[k][d],D)]
        for d = 1:D
            derivative = (bit[d] == 0 ? -1 : 1)/(hi[d]-lo[d])*_periodic_amr_ref_moment(Float64(ps.x[d][p]),Float64(clouds[d].half_width),lo[d],hi[d],Float64(mesh.L[d]),nothing)
            for j = 1:D
                j == d && continue
                derivative *= _periodic_amr_ref_moment(Float64(ps.x[j][p]),Float64(clouds[j].half_width),lo[j],hi[j],Float64(mesh.L[j]),bit[j])
            end
            field[d] -= raw[raw_index]*derivative
        end
    end
    return field
end
_periodic_amr_value(x) = x
_periodic_amr_value(x::Tuple) = map(_periodic_amr_value,x)
_periodic_amr_value(x::NamedTuple) = NamedTuple{keys(x)}(map(_periodic_amr_value,values(x)))
_periodic_amr_value(x::AbstractArray) = map(_periodic_amr_value,x)
function _periodic_amr_value(x::Union{AdaptiveCompositeESSolveReport,AdaptiveCompositeESStepReport,AdaptiveCompositeESRegridReport,HybridPlasmaPIC._AdaptiveCompositeCrossing,HybridPlasmaPIC._PeriodicCompositeTransport})
    return NamedTuple{fieldnames(typeof(x))}(ntuple(i -> _periodic_amr_value(getfield(x,i)),fieldcount(typeof(x))))
end
_periodic_amr_snapshot(st,species) = _periodic_amr_value(HybridPlasmaPIC._checkpoint_state(st,species))
function _periodic_amr_population(::Val{D},::Type{T};first_count = 1) where {D,T}
    a = ParticleSet{D,T}(first_count;q=T(1),m=T(2))
    b = ParticleSet{D,T}(1;q=T(-0.5),m=T(3))
    tracer = ParticleSet{D,T}(1;q=T(0.5),m=T(1))
    for p in eachindex(a.weight),d = 1:D
        a.x[d][p],a.v[d][p] = T((0.25,0.375,0.5)[d]),T((0.125,0.1875,0.0625)[d])
        a.weight[p] = T(1.25)/T(first_count)
    end
    for d = 1:D
        b.x[d][1],b.v[d][1] = T((2.75,3.25,2.5)[d]),T((-0.0625,0.03125,-0.09375)[d])
        tracer.x[d][1],tracer.v[d][1] = T((1.5,2.75,3.25)[d]),T(0.0625)
    end
    b.weight[1],tracer.weight[1] = T(2.5),zero(T)
    D == 2 && fill!(a.v[3],T(0.25))
    return [a,b,tracer,ParticleSet{D,T}(0)]
end

@testset "Periodic adaptive composite electrostatic essentials" begin
    H,Q = HybridPlasmaPIC,Rational{BigInt}
    for D in (2,3),T in (Float32,Float64)
        @testset "D=$D $T" begin
            base = CartesianAMRHierarchy(ntuple(_ -> 2,D),ntuple(_ -> T(4),D);periodic=true)
            fine = CartesianAMRHierarchy(base.nbase,base.L,((0,ntuple(_ -> 1,D)),);periodic=true)
            species = _periodic_amr_population(Val(D),T)
            clouds = ntuple(d -> d == 2 ? YeeTriangularCloud(Float32(0.5)) : YeeTriangularCloud(Float64(d == 1 ? 0.625 : 0.375)),D)
            st = AdaptiveCompositeElectrostaticPIC(fine,species;clouds,rtol=T(64)*eps(T))
            init!(st,species)
            reference = _periodic_amr_ref_field(fine,species,clouds,zero(Q))
            @test fine.periodic && fine.gauge == reference.gauge
            @test fine.masters == reference.masters && fine.free == reference.free
            @test isempty(fine.boundary) && !isempty(fine.seams)
            P = zeros(length(fine.vertices),length(fine.masters))
            for i in eachindex(fine.rows),(j,value) in fine.rows[i]
                P[i,j] = Float64(value)
            end
            @test isapprox(P,reference.P;rtol=1e-13,atol=1e-13)
            eigenvalues = eigvals(Symmetric(reference.K))
            @test abs(first(eigenvalues)) < 1e-12 && minimum(eigenvalues[2:end]) > 1e-5
            d = composite_electrostatic_diagnostics(st,species)
            @test iszero(sum(d.nodal_source)) && iszero(d.inventory.source_total)
            @test isapprox(Float64.(d.nodal_source),reference.source;rtol=2e-12,atol=2e-12)
            @test isapprox(Float64.(st.phi),reference.phi;rtol=4096eps(T),atol=4096eps(T))
            @test maximum(abs,d.full_gauss_residual) == d.solve_report.max_residual <= d.solve_report.tolerance
            @test d.gauge_residual == d.full_gauss_residual[fine.gauge] && all(iszero,d.harmonic_mean_field)
            u = zeros(length(fine.masters));u[fine.free] = Float64.(st.phi)
            @test isapprox(Float64(d.field_energy),dot(u,reference.K*u)/2;rtol=2e-12,atol=2e-12)
            forces = [_periodic_amr_ref_force(fine,reference,st.phi,ps,p,clouds) for ps in species for p in eachindex(ps.weight)]
            raw = H._acq_raw_potential(fine,st.phi)
            cursor = 0
            for ps in species,p in eachindex(ps.weight)
                cursor += 1
                actual = H._acc_gather(fine,raw,ntuple(j -> Q(ps.x[j][p]),D),clouds)
                @test isapprox(Float64.(collect(actual)),forces[cursor];rtol=2e-11,atol=2e-11)
            end
            @test composite_electric_field(st,ntuple(_ -> zero(T),D)) == composite_electric_field(st,base.L)
            saved0 = _periodic_amr_snapshot(st,species)
            @test saved0.integrator === :periodic_adaptive_composite_es
            @test keys(saved0.current) == (:refinements,:species,:phi,:last_solve_report,:last_step_report,:background_density,:transport)
            @test saved0.current.transport === nothing
            @test step!(st,species,-zero(T)) === st && regrid!(st,species,fine) === st
            @test isequal(_periodic_amr_snapshot(st,species),saved0)
            invalid = _periodic_amr_population(Val(D),T)
            invalid[1].weight[1] = one(T)
            @test_throws ArgumentError init!(st,invalid)
            @test isequal(_periodic_amr_snapshot(st,species),saved0)
            h = T(1/32)
            oldx,oldv = [map(copy,ps.x) for ps in species],[map(copy,ps.v) for ps in species]
            leaves = [Tuple(H._sce_particle_arrays(ps)) for ps in species]
            half = NTuple{3,T}[]
            cursor = 0
            for ps in species,p in eachindex(ps.weight)
                cursor += 1
                push!(half,ntuple(c -> T(Float64(ps.v[c][p])+Float64(h)/2*Float64(ps.q)/Float64(ps.m)*forces[cursor][c]),3))
            end
            step!(st,species,h)
            endpoint = _periodic_amr_ref_field(fine,species,clouds,zero(Q))
            @test isapprox(Float64.(st.phi),endpoint.phi;rtol=4096eps(T),atol=4096eps(T))
            cursor = 0
            for (s,ps) in enumerate(species),p in eachindex(ps.weight)
                cursor += 1
                end_force = _periodic_amr_ref_force(fine,endpoint,st.phi,ps,p,clouds)
                for axis = 1:D
                    expected = T(mod(Q(oldx[s][axis][p])+Q(h)*Q(half[cursor][axis]),Q(base.L[axis])))
                    @test isapprox(ps.x[axis][p],expected;rtol=64eps(T),atol=64eps(T))
                    @test ps.x[axis][p] != oldx[s][axis][p]
                end
                for c = 1:3
                    expected = T(Float64(half[cursor][c])+Float64(h)/2*Float64(ps.q)/Float64(ps.m)*end_force[c])
                    @test isapprox(ps.v[c][p],expected;rtol=64eps(T),atol=64eps(T))
                end
            end
            @test all(pair -> all(p -> p[1] === p[2],zip(pair[1],H._sce_particle_arrays(pair[2]))),zip(leaves,species))
            if D == 2
                @test species[1].v[3][1] == oldv[1][3][1]
            end
            report = st._state.step_report
            @test report.total_defect == report.field_defect+report.kick_defect+report.trajectory_defect
            @test !iszero(report.kick_work)
            @test st._state.transport !== nothing
            stepped = _periodic_amr_snapshot(st,species)
            @test_throws ArgumentError step!(st,species,-h)
            @test_throws ArgumentError regrid!(st,species[1:3],base)
            @test isequal(_periodic_amr_snapshot(st,species),stepped)
            regrid!(st,species,base)
            @test !iszero(st.last_regrid_report.mesh_work)
            history = _periodic_amr_value(st._history)
            regrid!(st,species,fine)
            @test isequal(_periodic_amr_value(st._history),history)
            @test st._history.state.refinements == fine.refinements && st._history.state.transport !== nothing
            @test iszero(composite_electrostatic_diagnostics(st,species).accounting_defect)
            original = H._checkpoint_state(st,species)
            saved = _periodic_amr_value(original)
            aliased = merge(original,(current=merge(original.current,(phi=st.phi,)),))
            @test_throws ArgumentError H._restore_checkpoint_state!(st,species,aliased)
            live_history = merge(original.history.state,(transport=st._history.state.transport,))
            @test_throws ArgumentError H._restore_checkpoint_state!(st,species,merge(original,(history=merge(original.history,(state=live_history,)),)))
            bad = deepcopy(original);bad.history.state.transport.defect[1][1][1] += 1
            @test_throws ArgumentError H._restore_checkpoint_state!(st,species,bad)
            @test_throws ArgumentError H._restore_checkpoint_state!(st,species,merge(original,(integrator=:adaptive_composite_es,)))
            @test isequal(_periodic_amr_snapshot(st,species),saved)
            mktempdir() do directory
                rawpath,archivepath = joinpath(directory,"periodic.chk"),joinpath(directory,"periodic.run")
                save_checkpoint(rawpath,st,species);archive_run(archivepath,st,species;rng_seed=0)
                other = _periodic_amr_population(Val(D),T;first_count=2)
                destination = AdaptiveCompositeElectrostaticPIC(base,other;clouds,rtol=T(64)*eps(T))
                init!(destination,other)
                for ps in other
                    foreach(x -> fill!(x,T(NaN)),ps.x)
                    fill!(ps.weight,T(NaN));ps.q,ps.m = T(NaN),-one(T)
                end
                fill!(destination.phi,T(NaN));fill!(destination._state.phi,T(NaN))
                destination.time,destination.step,destination.last_dt = T(Inf),-1,T(Inf)
                owner = other[1]
                load_checkpoint!(destination,other,rawpath)
                @test other[1] === owner && length(other[1].weight) == 1
                @test isequal(_periodic_amr_snapshot(destination,other),saved)
                step!(st,species,h);step!(destination,other,h)
                @test isequal(_periodic_amr_snapshot(destination,other),_periodic_amr_snapshot(st,species))
                load_archive!(destination,other,archivepath)
                @test isequal(_periodic_amr_snapshot(destination,other),saved)
                if D == 2 && T === Float64
                    adapt!(destination,other;refine_threshold=zero(T),coarsen_threshold=zero(T),max_level=2,buffer=1)
                    @test maximum(first,destination.mesh.addresses) == 2
                end
                adapt!(destination,other;refine_threshold=floatmax(T),coarsen_threshold=floatmax(T)/2,max_level=2,buffer=0)
                @test isempty(destination.mesh.refinements)
                @test iszero(composite_electrostatic_diagnostics(destination,other).accounting_defect)
            end
            # Nonneutral source is rejected unless its additional uniform charge
            # is explicitly requested and deposited, never removed from the RHS.
            one_species = [ParticleSet{D,T}(1)]
            for axis = 1:D
                one_species[1].x[axis][1] = T(0.125)
            end
            @test_throws ArgumentError AdaptiveCompositeElectrostaticPIC(base,one_species;clouds)
            wide = ntuple(_ -> YeeTriangularCloud(Float64(5.25)),D)
            neutralized = AdaptiveCompositeElectrostaticPIC(base,one_species;clouds=wide,background_charge=T(0.03125),periodic_neutrality=:neutralize,rtol=T(64)*eps(T))
            init!(neutralized,one_species)
            nd = composite_electrostatic_diagnostics(neutralized,one_species)
            rho = -one(Q)/prod(Q,base.L)
            @test nd.background_density == rho && nd.inventory.declared_background == Q(T(0.03125))*prod(Q,base.L)
            @test nd.inventory.neutralizing_charge == -1-nd.inventory.declared_background && nd.inventory.background == -1 && iszero(nd.inventory.source_total)
            wide_ref = _periodic_amr_ref_field(base,one_species,wide,rho)
            @test isapprox(Float64.(nd.nodal_source),wide_ref.source;rtol=2e-12,atol=2e-12)
            if D == 2 && T === Float64
                uniform_clouds = ntuple(d -> YeeTriangularCloud(T(8)*base.L[d]),D)
                uniform = AdaptiveCompositeElectrostaticPIC(base,one_species;clouds=uniform_clouds,periodic_neutrality=:neutralize)
                init!(uniform,one_species)
                @test all(iszero,composite_electrostatic_diagnostics(uniform,one_species).nodal_source)
            end
        end
    end
    # Same wrapped endpoints do not erase positive/negative multi-period paths.
    for D in (2,3),T in (Float32,Float64)
        mesh = CartesianAMRHierarchy(ntuple(_ -> 2,D),ntuple(_ -> T(4),D);periodic=true)
        neutral = ParticleSet{D,T}(1;q=zero(T))
        for axis = 1:D
            neutral.x[axis][1] = axis == 2 ? T(3.5) : T(0.5)
            neutral.v[axis][1] = axis == 2 ? T(-8) : T(8)
        end
        clouds = ntuple(_ -> YeeTriangularCloud(T(0.375)),D)
        st = AdaptiveCompositeElectrostaticPIC(mesh,[neutral];clouds)
        init!(st,[neutral]);step!(st,[neutral],one(T))
        d = composite_electrostatic_diagnostics(st,[neutral])
        @test length(d.last_step_report.crossings) == 4
        @test all(event -> all(event.axes),d.last_step_report.crossings)
        @test [event.parameter for event in d.last_step_report.crossings] == [3//16,7//16,11//16,15//16]
        for axis = 1:D
            @test d.transport.winding[1][axis][1] == (axis == 2 ? -2 : 2)
            @test neutral.x[axis][1] == (axis == 2 ? T(3.5) : T(0.5))
            @test iszero(d.transport.defect[1][axis][1])
        end
        @test count(i -> d.transport.tiles_before[i] != d.transport.tiles_after[i],eachindex(d.transport.tiles_before)) == 2
        # Ideal remains below L, but its nearest-T endpoint is L and must fold to0.
        for axis = 1:D
            neutral.x[axis][1] = prevfloat(mesh.L[axis])
            neutral.v[axis][1] = mesh.L[axis]-prevfloat(mesh.L[axis])
        end
        init!(st,[neutral]);step!(st,[neutral],T(0.5))
        d = composite_electrostatic_diagnostics(st,[neutral])
        @test only(d.last_step_report.crossings).parameter == 1
        for axis = 1:D
            @test iszero(neutral.x[axis][1]) && d.transport.winding[1][axis][1] == 1
            @test d.transport.defect[1][axis][1] == (Q(mesh.L[axis])-Q(prevfloat(mesh.L[axis])))/2
        end
    end
    huge = H._periodic_represented_lift((BigInt(1)<<70)*Q(4)+1//8,4.0)
    @test huge.wrapped == 0.125 && huge.winding == (BigInt(1)<<70) && iszero(huge.defect)

    @testset "Full residual, real solver failure and preserved ground schema" begin
        T = Float64
        mesh = CartesianAMRHierarchy((2,2),(T(4),T(4));periodic=true)
        clouds = (YeeTriangularCloud(T(0.375)),YeeTriangularCloud(T(0.375)))
        particles = [ParticleSet{2,T}(1),ParticleSet{2,T}(1;q=-one(T))]
        particles[1].x[1][1],particles[1].x[2][1] = T(0.5),T(1)
        particles[2].x[1][1],particles[2].x[2][1] = T(1),T(0.5)
        st = AdaptiveCompositeElectrostaticPIC(mesh,particles;clouds,rtol=T(1024)*eps(T))
        init!(st,particles)
        d = composite_electrostatic_diagnostics(st,particles)
        @test iszero(d.nodal_source[mesh.gauge])
        reference = _periodic_amr_ref_field(mesh,particles,clouds,zero(Q))
        target = d.solve_report.tolerance
        badphi = T.(reference.K[mesh.free,mesh.free]\(Float64.(d.nodal_source[mesh.free]).+Float64(3target/4)))
        residual,free = H._acq_residual(st._state.operator,badphi,d.nodal_source)
        free_max,full_max = maximum(abs,free),maximum(abs,residual)
        @test free_max < target < full_max
        saved = H._checkpoint_state(st,particles)
        receipt = d.solve_report
        wrong_report = AdaptiveCompositeESSolveReport{T}(receipt.iterations,receipt.status,true,free_max,receipt.source_norm,receipt.tolerance,receipt.scale)
        candidate = merge(saved,(current=merge(saved.current,(phi=badphi,last_solve_report=wrong_report)),initial_energy=H._kinetic_energy_rational(particles)+H._acq_energy(st._state.operator,badphi)))
        before = _periodic_amr_value(saved)
        @test_throws ArgumentError H._restore_checkpoint_state!(st,particles,candidate)
        @test isequal(_periodic_amr_snapshot(st,particles),before)

        # Exact zero-source initialization succeeds with maxiter1. The moved
        # source has reduced RHS[+.25,-.25,-.25], which is not a CG eigenvector.
        coincident = [ParticleSet{2,T}(1),ParticleSet{2,T}(1;q=-one(T))]
        moved = [ParticleSet{2,T}(1),ParticleSet{2,T}(1;q=-one(T))]
        for ps in coincident
            ps.x[1][1],ps.x[2][1] = T(0.5),T(1)
        end
        moved[1].x[1][1],moved[1].x[2][1] = T(0.5),T(1)
        moved[2].x[1][1],moved[2].x[2][1] = T(2.5),T(3)
        limited = AdaptiveCompositeElectrostaticPIC(mesh,coincident;clouds,maxiter=1,rtol=T(64)*eps(T))
        init!(limited,coincident)
        before = _periodic_amr_snapshot(limited,coincident)
        @test_throws ArgumentError init!(limited,moved)
        @test isequal(_periodic_amr_snapshot(limited,coincident),before)
        step!(limited,coincident,T(1/8))
        @test limited.step == 1

        # The old scalar Yee cloud supplies an independent agreement surface
        # on its exact-tiling, narrow stationary-image admitted domain.
        one_species = [ParticleSet{2,T}(1)]
        one_species[1].x[1][1],one_species[1].x[2][1] = T(0.125),T(0.125)
        agreement = AdaptiveCompositeElectrostaticPIC(mesh,one_species;clouds,periodic_neutrality=:neutralize)
        grid = YeeGrid((2,),(T(4),);boundaries=(:periodic,))
        for position in (T(0.125),T(3.875))
            one_species[1].x[1][1],one_species[1].x[2][1] = position,position
            init!(agreement,one_species)
            coefficients = fill(zero(Q),2)
            H._pyee_each_basis(Q(position),Q(position),grid,Q(0.375),Val(true)) do index,value
                coefficients[index] += value
            end
            expected = [prod(axis -> coefficients[Int(mesh.vertices[raw][axis]/2)+1],1:2)-1//4 for raw in mesh.masters]
            @test composite_electrostatic_diagnostics(agreement,one_species).nodal_source == expected
        end

        groundmesh = CartesianAMRHierarchy((2,2),(T(4),T(4)))
        ground_particles = [ParticleSet{2,T}(1)]
        ground_particles[1].x[1][1],ground_particles[1].x[2][1] = T(1.25),T(1.5)
        @test_throws ArgumentError AdaptiveCompositeElectrostaticPIC(groundmesh,ground_particles;clouds,periodic_neutrality=:neutralize)
        ground = AdaptiveCompositeElectrostaticPIC(groundmesh,ground_particles;clouds)
        init!(ground,ground_particles)
        legacy = H._checkpoint_state(ground,ground_particles)
        @test legacy.integrator === :adaptive_composite_es
        @test keys(legacy.current) == (:refinements,:species,:phi,:last_solve_report,:last_step_report)
        @test keys(legacy.configuration) == (:D,:T,:nbase,:lengths,:ratio,:ghosts,:clouds,:background,:factors,:rtol,:atol,:maxiter,:boundary,:field,:solve,:force,:phase,:adaptation,:arithmetic,:events)
        @test legacy.configuration.boundary === :grounded_closest_point_bound_tails_v1
        @test legacy.configuration.phase === :synchronous_global_fixed_duration_v1
        step!(ground,ground_particles,T(1/64))
        @test fieldnames(typeof(ground._state.step_report)) == (:dt,:time_before,:time_after,:step_before,:step_after,:epoch_start_time,:epoch_start_step,:refinements,:kinetic_before,:kinetic_after,:field_before,:field_after,:source_work,:kick_work,:field_defect,:kick_defect,:trajectory_defect,:total_defect,:max_residual,:crossings)
        @test fieldnames(typeof(ground._state.solve_report)) == (:iterations,:status,:solved,:max_residual,:source_norm,:tolerance,:scale)
    end
end
