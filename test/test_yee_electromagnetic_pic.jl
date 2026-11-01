module TestYeeElectromagneticPIC
using HybridPlasmaPIC, Test, LinearAlgebra, Serialization
import JET
include("yee_charge_oracle.jl")
include("sourced_yee_oracle.jl")
const O = YeeChargeOracle
const F = SourcedYeeOracle
const H = HybridPlasmaPIC
const Q = Rational{BigInt}

function fixture(::Type{T}; collection = false, transverse = true, width = T(0.19)) where {T}
    g = YeeGrid((8,), (one(T),))
    ps = ParticleSet{1,T}(2; q = T(0.3), m = T(1.2))
    ps.x[1] .= T[0.34, 0.63]
    ps.weight .= T[0.02, 0.03]
    ps.v[1] .= T[0.012, -0.017]
    if transverse
        ps.v[2] .= T[0.04, -0.02]
        ps.v[3] .= T[-0.03, 0.025]
    end
    st = YeeElectromagneticPIC(g; shape = YeeTriangularCloud(width), c = T(0.7),
        background_charge = T(-0.015), left_plate_charge = T(0.01))
    if transverse
        for i = 1:9
            x = T(i - 1) / 8
            st.fields.E[2][i] = T(0.04) * sinpi(x)
            st.fields.E[3][i] = T(0.03) * sinpi(2x)
        end
        st.fields.E[2][end] = st.fields.E[3][end] = zero(T)
        st.fields.B[2] .= T(0.025)
        st.fields.B[3] .= T(-0.02)
    end
    species = collection ? [ps, ParticleSet{1,T}(0; q = T(-0.2), m = T(2))] : ps
    init!(st, species)
    return st, species
end
pops(ps::ParticleSet) = (ps,)
pops(species::AbstractVector) = species
function snapshot(st, particles)
    return deepcopy(([(p.x,p.v,p.weight,p.id,p.tag,p.q,p.m) for p in pops(particles)],
        st.fields.E,st.fields.B,st.time[],st.step[],st.initialized[],st.last_dt[],
        st.right_plate_charge[],st.last_impulse,st.last_work[],
        st.last_modified_energy_defect[],st.last_physical_energy_change[]))
end

function energy_reference(st,particles,fields)
    kinetic = sum((Q(p.m)*Q(p.weight[i])*sum(Q(p.v[c][i])^2 for c=1:3)/2
        for p in pops(particles) for i in eachindex(p.weight));init=zero(Q))
    electric = sum(Q(x)^2 for a in fields.E for x in a)
    magnetic = sum(Q(x)^2 for a in fields.B for x in a)
    return kinetic+Q(st.g.dx[1])*(electric+Q(st.c)^2*magnetic)/2
end
replace_field(st,name,value) = typeof(st)((field===name ? value : getfield(st,field)
    for field in fieldnames(typeof(st)))...)

# Independent physical-volume/path quadrature. No production A/g/U/rounding
# helper contributes to this oracle. The existing quadrature is exact for the
# piecewise polynomials, split at physical nodes and cloud breakpoints.
function pathnodes(x, y, g, a)
    x, y, a, L = Q(x), Q(y), Q(a), Q(g.L[1])
    x == y && return O.cloud_nodes(x,g.n[1],L,a)
    lo, hi = minmax(x,y)
    knots = Q[lo,hi]
    for node = 0:g.n[1], shift in (-a,zero(Q),a)
        knot = node * L / g.n[1] + shift
        lo < knot < hi && push!(knots,knot)
    end
    return [O.polynomial_integral(s -> O.cloud_nodes(s,g.n[1],L,a)[j],knots)/(hi-lo)
        for j = 1:g.n[1]+1]
end
function gather(E,B,x,st)
    n,L,a = st.g.n[1],st.g.L[1],st.shape.half_width
    C = O.cloud_nodes(x,n,L,a)
    edge(field) = Q(st.g.dx[1]) * sum(Q(field[e+1])*O.cloud_edge(x,e,n,L,a) for e = 0:n-1)
    return Q[edge(E[1]),sum(Q.(E[2]).*C),sum(Q.(E[3]).*C)], Q[0,edge(B[2]),edge(B[3])]
end
# Independent implicit-midpoint matrix representation of frozen-field Boris.
function kick(v,E,B,q,m,duration)
    K = Q[0 B[3] -B[2]; -B[3] 0 B[1]; B[2] -B[1] 0]
    a = Q(q)/Q(m)*Q(duration)/2
    eye = Matrix{Q}(I,3,3)
    return (eye-a*K) \ ((eye+a*K)*Q.(v)+2 * a * E)
end
function reference(st,particles,dt)
    T = eltype(st.g.L)
    species = pops(particles)
    u = [ntuple(c -> zeros(T,length(ps.weight)),3) for ps in species]
    y = [zeros(T,length(ps.weight)) for ps in species]
    for (k,ps) in enumerate(species), p in eachindex(ps.weight)
        E,B = gather(st.fields.E,st.fields.B,ps.x[1][p],st)
        velocity = kick([ps.v[c][p] for c=1:3],E,B,ps.q,ps.m,Q(dt)/2)
        for c=1:3
            u[k][c][p] = T(velocity[c])
        end
        y[k][p] = T(Q(ps.x[1][p])+Q(dt)*Q(u[k][1][p]))
    end
    current = ntuple(c -> zeros(Q,length(st.fields.E[c])),3)
    current[1] .= O.cloud_impulse(species,[ps.x[1] for ps in species],y,
        st.g.n[1],st.g.L[1],st.shape.half_width)
    for (k,ps) in enumerate(species), p in eachindex(ps.weight)
        C = pathnodes(ps.x[1][p],y[k][p],st.g,st.shape.half_width)
        for c=2:3, j=2:st.g.n[1]
            current[c][j] += Q(ps.q)*Q(ps.weight[p])*Q(dt)*Q(u[k][c][p])*C[j]/Q(st.g.dx[1])
        end
    end
    impulse = map(a -> T.(a),current)
    speed = st isa YeeElectromagneticPIC ? st.c : one(T)
    En,Bn,_ = F.step(st.fields.E,st.fields.B,impulse,st.g.n[1],false,st.g.dx[1],dt,speed)
    v = [ntuple(c -> zeros(T,length(ps.weight)),3) for ps in species]
    for (k,ps) in enumerate(species), p in eachindex(ps.weight)
        E,B = gather(En,Bn,y[k][p],st)
        velocity = kick([u[k][c][p] for c=1:3],E,B,ps.q,ps.m,Q(dt)/2)
        for c=1:3
            v[k][c][p] = T(velocity[c])
        end
    end
    return y,v,impulse,En,Bn
end

@testset "Independent nodal path and PEC source pairing" begin
    for T in (Float32,Float64), W in (Float32,Float64)
        st,ps = fixture(T;width=W(0.19))
        h,a = Q(st.g.L[1])/8,Q(st.shape.half_width)
        for (x,y) in ((T(0.01),T(0.71)),(T(0.71),T(0.01)),(T(0.5),T(0.5)))
            expected = pathnodes(x,y,st.g,a)
            for node=1:7
                @test H._yee_cloud_node_path(Q(x),Q(y),node,h,a) == expected[node+1]
            end
        end
        for x in T[0.001,0.125,0.5,0.999]
            E,B = gather(st.fields.E,st.fields.B,x,st)
            actual = H._ye_gather(st.fields,x,st.g,st.shape)
            @test actual[1] == Tuple(E)
            @test actual[2] == Tuple(B)
        end
        candidate = H._yl_make_stage([ps],[ps],T).candidates
        candidate[1].x[1] .= T[0.01,0.81]
        current = map(a -> zero.(a),st.last_impulse)
        H._yee_transverse_impulse!(current,[ps],candidate,st.g,st.shape,T(0.01))
        for c=2:3
            @test iszero(current[c][1]) && iszero(current[c][end])
            for j=2:8
                expected = sum(Q(ps.q)*Q(ps.weight[p])*Q(T(0.01))*Q(candidate[1].v[c][p])*
                    pathnodes(ps.x[1][p],candidate[1].x[1][p],st.g,a)[j]/Q(st.g.dx[1]) for p=1:2)
                @test O.rounded(current[c][j],expected)
            end
        end
        # Whole-entry signed cancellation retains sources whose q*w would Inf.
        sources = [ParticleSet{1,T}(1;q=q) for q in (floatmax(T),-floatmax(T),one(T))]
        for (i,p) in enumerate(sources)
            p.x[1] .= T(0.5); p.weight .= i<3 ? T(2) : one(T)
            p.v[2] .= one(T); p.v[3] .= -one(T)
        end
        candidates = H._yl_make_stage(sources,sources,T).candidates
        H._yee_transverse_impulse!(current,sources,candidates,st.g,st.shape,T(0.01))
        C = pathnodes(T(0.5),T(0.5),st.g,a)
        for j=2:8
            expected = Q(T(0.01))*C[j]/Q(st.g.dx[1])
            @test O.rounded(current[2][j],expected)
            @test O.rounded(current[3][j],-expected)
        end
    end
end

@testset "Full 1D3V map, continuity and physical/modified diagnostics" begin
    for T in (Float32,Float64)
        st,ps = fixture(T;collection=true)
        dt = T(0.005)
        expected = reference(st,ps,dt)
        before = snapshot(st,ps)
        identities = map(objectid,(st.fields.E...,st.fields.B...,ps[1].x...,ps[1].v...))
        oldK = kinetic_energy(ps)
        oldH = energy_reference(st,ps,st.fields)
        oldmodified = yee_modified_energy(st,dt)
        oldcharge = O.cloud_charge(ps,[copy(p.x[1]) for p in ps],8,st.g.L[1],st.shape.half_width)
        step!(st,ps,dt)
        for (k,p) in enumerate(ps)
            @test p.x[1] ≈ expected[1][k] rtol=64eps(T) atol=64eps(T)
            for c=1:3
                @test p.v[c] ≈ expected[2][k][c] rtol=128eps(T) atol=128eps(T)
            end
        end
        for c=1:3
            @test st.last_impulse[c] ≈ expected[3][c] rtol=128eps(T) atol=128eps(T)
            @test st.fields.E[c] ≈ expected[4][c] rtol=128eps(T) atol=128eps(T)
            @test st.fields.B[c] ≈ expected[5][c] rtol=128eps(T) atol=128eps(T)
        end
        newcharge = O.cloud_charge(ps,[p.x[1] for p in ps],8,st.g.L[1],st.shape.half_width)
        @test maximum(abs,newcharge-oldcharge+O.incidence(st.last_impulse[1])) <= 64 * Q(eps(T))
        @test identities == map(objectid,(st.fields.E...,st.fields.B...,ps[1].x...,ps[1].v...))
        d = yee_electromagnetic_diagnostics(st,ps)
        @test O.rounded(d.last_modified_energy_defect,Q(kinetic_energy(ps))-Q(oldK)-Q(st.last_work[]))
        @test O.rounded(d.last_physical_energy_change,energy_reference(st,ps,st.fields)-oldH)
        @test abs(Q(yee_modified_energy(st,dt))-Q(oldmodified)+Q(st.last_work[])) <= 256 * Q(eps(T))
        @test any(!iszero,st.last_impulse[2]) && any(!iszero,st.last_impulse[3])
        @test !isequal(snapshot(st,ps),before)
    end
end

@testset "Longitudinal reduction, original half duration and temporal refinement" begin
    for T in (Float32,Float64)
        source,_ = fixture(T)
        empty = ParticleSet{1,T}[]
        vacuumPIC = YeeElectromagneticPIC(source.g;shape=source.shape,c=source.c)
        foreach(copyto!,vacuumPIC.fields.E,source.fields.E)
        foreach(copyto!,vacuumPIC.fields.B,source.fields.B)
        init!(vacuumPIC,empty)
        vacuum = MaxwellVacuum(source.g;c=source.c)
        foreach(copyto!,vacuum.fields.E,vacuumPIC.fields.E)
        foreach(copyto!,vacuum.fields.B,vacuumPIC.fields.B)
        step!(vacuumPIC,empty,T(0.005)); step_maxwell!(vacuum,T(0.005))
        @test isequal(vacuumPIC.fields.E,vacuum.fields.E)
        @test isequal(vacuumPIC.fields.B,vacuum.fields.B)
        st,ps = fixture(T;transverse=false)
        longitudinal = YeeLongitudinalPIC(st.g;shape=st.shape,background_charge=st.background_charge,
            left_plate_charge=st.left_plate_charge)
        other = deepcopy(ps)
        init!(longitudinal,other)
        for _=1:2
            step!(st,ps,T(0.005)); step!(longitudinal,other,T(0.005))
            @test isequal(st.fields.E,longitudinal.fields.E)
            @test isequal(st.fields.B,longitudinal.fields.B)
            @test isequal(ps.x,other.x) && isequal(ps.v,other.v)
        end
        # Odd subnormal half duration, with an unrepresentable q/m rescued by
        # its original field/duration factors. Tracers do not perturb fields.
        m = floatmin(T)
        tracer = ParticleSet{1,T}(1;q=ldexp(one(T),-exponent(m)),m=m)
        tracer.x[1] .= T(0.25); tracer.weight .= zero(T)
        thin = YeeElectromagneticPIC(YeeGrid((2,),(one(T),));shape=YeeTriangularCloud(T(0.125)),left_plate_charge=m)
        init!(thin,tracer)
        tau = nextfloat(zero(T))
        step!(thin,tracer,tau)
        @test O.rounded(tracer.v[1][1],Q(tau)*Q(tracer.q))
    end
    # Coupled, nonzero transverse evolution: compare common integer times at
    # successive refinements. No assumed exactly conserved discrete energy.
    function final(h)
        st,ps = fixture(Float64)
        initial = energy_reference(st,ps,st.fields)
        for _=1:round(Int,0.04/h)
            step!(st,ps,h)
        end
        values = vcat(ps.x[1],ps.v...,st.fields.E...,st.fields.B...)
        return values,Float64(energy_reference(st,ps,st.fields)-initial)
    end
    coarse,medium,fine = final(0.004),final(0.002),final(0.001)
    @test 3 < norm(coarse[1]-medium[1])/norm(medium[1]-fine[1]) < 5
    @test abs(fine[2]) < abs(medium[2]) < abs(coarse[2])
end

@testset "Actual coupled transverse rest map" begin
    for T in (Float32,Float64)
        g = YeeGrid((4,),(one(T),))
        st = YeeElectromagneticPIC(g;shape=YeeTriangularCloud(T(0.17)),c=T(0.7))
        particles = [ParticleSet{1,T}(1;q=q) for q in (one(T),-one(T))]
        for p in particles
            p.x[1] .= T(0.5); p.weight .= T(0.02)
        end
        amplitude = T === Float32 ? T(0.001) : T(0.00001)
        st.fields.E[2][2:4] .= amplitude .* T[1,0.5,-0.3]
        init!(st,particles)
        delta = Q(g.dx[1])
        D = zeros(Q,4,3)
        for edge=1:4,node=1:3
            D[edge,node] = ((node==edge ? 1 : 0)-(node==edge-1 ? 1 : 0))/delta
        end
        C = O.cloud_nodes(T(0.5),4,one(T),st.shape.half_width)[2:4]
        plasma = sum(Q(p.q)^2*Q(p.weight[1])/Q(p.m) for p in particles)/delta
        A = Q(st.c)^2*(D'*D)+plasma*(C*C')
        dt = T(0.02)
        initial = Q.(st.fields.E[2][2:4])
        expected = initial-Q(dt)^2/2 .* (A*initial)
        @test Q(dt)^2*maximum(eigvals(Symmetric(Float64.(A)))) < 4
        step!(st,particles,dt)
        @test st.fields.E[2][2:4] ≈ T.(expected) rtol=64eps(T) atol=64eps(T)*amplitude
        @test any(!iszero,st.fields.B[3]) && any(!iszero,particles[1].v[2])
    end
end

@testset "Original geometry through both coupled electric maps" begin
    T = Float32
    eta = nextfloat(zero(T))
    # Zero-weight tracers isolate the force lifecycle: fields/current remain
    # unchanged, but both rounded velocity stages and the stored drift matter.
    cases = (
        (name=:underflow,n=2,a=ldexp(1.0,100),field=eta,x=T(0.25),v=zero(T),
            q=ldexp(one(T),126),m=ldexp(one(T),-126),dt=ldexp(one(T),-56)),
        (name=:overflow,n=7,a=ldexp(1.0,-10),field=floatmax(T),x=T(0.5),v=zero(T),
            q=ldexp(one(T),-126),m=ldexp(one(T),126),dt=ldexp(one(T),-20)),
        (name=:cancellation,n=7,a=ldexp(1.0,-10),field=one(T),x=T(0.5),v=-one(T),
            q=T(256),m=one(T),dt=T(1/128)),
    )
    for Model in (YeeLongitudinalPIC,YeeElectromagneticPIC), case in cases
        g = YeeGrid((case.n,),(one(T),))
        ps = ParticleSet{1,T}(1;q=case.q,m=case.m)
        ps.x[1][1],ps.v[1][1],ps.weight[1] = case.x,case.v,zero(T)
        st = Model(g;shape=YeeTriangularCloud(case.a),left_plate_charge=case.field)
        init!(st,ps)
        E,B = gather(st.fields.E,st.fields.B,case.x,st)
        first = kick([case.v,zero(T),zero(T)],E,B,ps.q,ps.m,Q(case.dt)/2)
        if case.name === :underflow
            @test 0 < E[1] < Q(eta)/2
            @test yee_longitudinal_field(st.fields.E[1],case.x,g;shape=st.shape) === zero(T)
            @test Q(ps.q)/Q(ps.m)*Q(case.dt)*E[1] > Q(floatmin(T))
        elseif case.name === :overflow
            midpoint = Q(floatmax(T))+Q(2)^(exponent(floatmax(T))-precision(T))
            @test Q(7)*Q(g.dx[1]) == 1+3Q(2)^(-26)
            @test E[1]-midpoint == Q(2)^78*(Q(2)^24-3) > 0
            @test ps.q/ps.m === zero(T)
            @test_throws ArgumentError yee_longitudinal_field(st.fields.E[1],case.x,g;shape=st.shape)
            @test 0 < first[1] < Q(floatmin(T))
        else
            @test E[1] == 1+3Q(2)^(-26)
            @test first[1] == 3Q(2)^(-26)
        end
        expected = reference(st,ps,case.dt)
        step!(st,ps,case.dt)
        @test isequal(ps.x[1],expected[1][1])
        @test isequal(ps.v,expected[2][1])
        @test isequal(st.last_impulse,expected[3])
        @test isequal(st.fields.E,expected[4]) && isequal(st.fields.B,expected[5])
        @test ps.x[1][1] === case.x && !iszero(ps.v[1][1])
        if case.name === :cancellation
            @test ps.v[1][1] === nextfloat(one(T))
        end
    end
end

@testset "Unrepresentable magnetic gather has a finite Cayley rotation" begin
    T = Float32
    g = YeeGrid((7,),(one(T),))
    ps = ParticleSet{1,T}(1;q=ldexp(one(T),-100))
    ps.x[1][1],ps.weight[1] = T(0.5),zero(T)
    ps.v[1][1],ps.v[3][1] = one(T),one(T)
    st = YeeElectromagneticPIC(g;shape=YeeTriangularCloud(ldexp(1.0,-10)))
    st.fields.B[2] .= floatmax(T)
    init!(st,ps)
    E,B = gather(st.fields.E,st.fields.B,ps.x[1][1],st)
    midpoint = Q(floatmax(T))+Q(2)^(exponent(floatmax(T))-precision(T))
    @test B[2] > midpoint && all(iszero,E)
    dt = ldexp(one(T),-40)
    expected = reference(st,ps,dt)
    before = copy(ps.v[3])
    step!(st,ps,dt)
    @test isequal(ps.v,expected[2][1]) && ps.v[3] != before
    @test isequal(st.fields.E,expected[4]) && isequal(st.fields.B,expected[5])
end

@testset "Final rational half-kick overflow rolls back and permits retry" begin
    for T in (Float32,Float64), Model in (YeeLongitudinalPIC,YeeElectromagneticPIC)
        g = YeeGrid((2,),(one(T),))
        ps = ParticleSet{1,T}(1;q=ldexp(one(T),-exponent(floatmin(T))),m=floatmin(T))
        ps.x[1][1],ps.v[1][1],ps.weight[1] = T(0.5),prevfloat(floatmax(T)),zero(T)
        st = Model(g;shape=YeeTriangularCloud(T(1/8)),left_plate_charge=T(4))
        init!(st,ps)
        dt = nextfloat(zero(T))
        alpha = Q(ps.q)/Q(ps.m)
        ulp = Q(floatmax(T))-Q(prevfloat(floatmax(T)))
        @test alpha*Q(dt)*4/2 == ulp
        # For this uniform N2 state, the EM dyadic field bound is exactly 8;
        # it also bounds the longitudinal force. Stiffness is zero.
        @test Q(dt)*Q(ps.v[1][1])+Q(dt)^2*alpha*8/2 < 1//2
        @test Q(dt) < Q(g.dx[1])
        E,B = gather(st.fields.E,st.fields.B,ps.x[1][1],st)
        first = T.(kick([ps.v[c][1] for c=1:3],E,B,ps.q,ps.m,Q(dt)/2))
        @test first[1] === floatmax(T)
        y = T(Q(ps.x[1][1])+Q(dt)*Q(first[1]))
        @test Q(1)/8 < Q(y) < Q(7)/8
        En,Bn = gather(st.fields.E,st.fields.B,y,st)
        final = kick(first,En,Bn,ps.q,ps.m,Q(dt)/2)
        @test final[1] == Q(2)^(exponent(floatmax(T))+1)
        before = deepcopy(H._checkpoint_state(st,ps))
        @test_throws ArgumentError step!(st,ps,dt)
        @test isequal(H._checkpoint_state(st,ps),before)
        ps.v[1][1] = zero(T)
        expected = reference(st,ps,dt)
        step!(st,ps,dt)
        @test isequal(ps.v,expected[2][1]) && st.step[] == 1
    end
end

@testset "Admission and late failure preserve public state" begin
    for T in (Float32,Float64)
        st,ps = fixture(T)
        before = snapshot(st,ps)
        @test_throws ArgumentError step!(st,ps,maxwell_cfl(st))
        @test isequal(snapshot(st,ps),before)
        aliasfields = YeeFields{1,T}((st.fields.E[1],st._candidate.E[2],st._candidate.E[3]),st._candidate.B)
        aliased = replace_field(st,:_candidate,aliasfields)
        @test_throws ArgumentError step!(aliased,ps,T(0.001))
        @test isequal(snapshot(st,ps),before)
        # Source mutation requires explicit reinitialization.
        ps.weight[1] = nextfloat(ps.weight[1])
        invalid = snapshot(st,ps)
        @test_throws ArgumentError step!(st,ps,T(0.001))
        @test isequal(snapshot(st,ps),invalid)
        ps.weight[1] = before[1][1][3][1]
        # Combined rest boundary: independent vacuum and plasma checks both
        # pass, while their composed sum fails before publication.
        tau = T(0.9)*maxwell_cfl(st)
        old = copy(ps.weight)
        ps.weight .*= T(1000)
        init!(st,ps)
        energy = H._ye_physical_energy(st,(ps,),st.fields)
        plasma = sum(Q(ps.q)^2*Q(w)/Q(ps.m) for w in ps.weight)/Q(st.g.dx[1])
        vacuum = 4 * Q(st.c)^2/Q(st.g.dx[1])^2
        @test Q(tau)^2*plasma < 4 && Q(tau)^2*vacuum < 4
        @test Q(tau)^2*(plasma+vacuum) > 4
        @test_throws ArgumentError H._ye_admit_step(st,(ps,),tau,energy)
        ps.weight .= old; init!(st,ps)
        # Full-speed tube catches transverse motion although vx=0; this is a
        # physical exclusion guard, independently of the stored straight path.
        ps.v[1] .= zero(T); ps.v[2] .= T(100)
        init!(st,ps)
        before = snapshot(st,ps)
        @test_throws ArgumentError step!(st,ps,T(0.01))
        @test isequal(snapshot(st,ps),before)
        # An invalid last diagnostic is reached after staging even for q=0;
        # physical energy itself need not round finite for tube construction.
        large = T === Float32 ? 80 : 600
        tracer = ParticleSet{1,T}(1;q=zero(T))
        tracer.x[1] .= T(0.5); tracer.v[2] .= ldexp(one(T),large)
        late = YeeElectromagneticPIC(st.g;shape=st.shape)
        init!(late,tracer)
        before = snapshot(late,tracer)
        @test_throws ArgumentError step!(late,tracer,ldexp(one(T),-large-20))
        @test isequal(snapshot(late,tracer),before)
        for (mode,flush) in ((RoundDown,false),(RoundNearest,true))
            saved = Base.Rounding.rounding_raw(Float64),get_zero_subnormals()
            before = snapshot(st,ps)
            try
                Base.Rounding.setrounding_raw(Float64,Base.Rounding.to_fenv(mode))
                Base.set_zero_subnormals(flush)
                @test_throws ArgumentError step!(st,ps,T(0.001))
            finally
                Base.Rounding.setrounding_raw(Float64,saved[1]); Base.set_zero_subnormals(saved[2])
            end
            @test isequal(snapshot(st,ps),before)
        end
    end
end

@testset "Electromagnetic checkpoint/archive and concrete compiler contracts" begin
    for T in (Float32,Float64), collection in (false,true)
        st,ps = fixture(T;collection)
        step!(st,ps,T(0.002))
        # Preserve actual history that rebuilding/projecting would erase.
        st.fields.E[1][2] = nextfloat(st.fields.E[1][2])
        st.fields.B[2][2] = nextfloat(st.fields.B[2][2])
        mktempdir() do dir
            path,archive = joinpath(dir,"em.bin"),joinpath(dir,"archive.bin")
            save_checkpoint(path,st,ps); archive_run(archive,st,ps;rng_seed=23)
            saved = deserialize(path)
            @test saved.integrator === :yee_electromagnetic_pic && saved.schema_version === 2
            @test saved.force_policy === :original_geometry_boris_nearest_even
            @test saved.force_policy_version === 1
            @test isequal(saved,load_archive(archive).state)
            restored,_ = fixture(T;collection)
            viewof(::Type{A},n) where {A} = view(zeros(A,2n),1:2:2n)
            viewed = ParticleSet{1,T}((viewof(T,2),),ntuple(_ -> viewof(T,2),3),
                viewof(T,2),viewof(UInt64,2),viewof(UInt32,2),zero(T),one(T))
            other = collection ? [viewed,ParticleSet{1,T}(0)] : viewed
            particleids = H._yl_particle_arrays(viewed)
            ids = map(objectid,(restored.fields.E...,restored.fields.B...))
            load_checkpoint!(restored,other,path)
            @test isequal(snapshot(restored,other),snapshot(st,ps))
            @test ids == map(objectid,(restored.fields.E...,restored.fields.B...))
            @test all(pair -> pair[1]===pair[2],zip(particleids,H._yl_particle_arrays(viewed)))
            @test all(candidate -> all(a -> a isa Vector, H._yl_particle_arrays(candidate)),
                restored._stage[].candidates)
            @test all(pair -> !Base.mightalias(pair...),
                zip(H._yl_particle_arrays(restored._stage[].candidates[1]),H._yl_particle_arrays(viewed)))
            # Empty ordinary particle storage may resize without replacing arrays.
            ordinary = ParticleSet{1,T}(0)
            destination = collection ? [ordinary,ParticleSet{1,T}(0)] : ordinary
            fresh,_ = fixture(T;collection)
            ordinaryids = H._yl_particle_arrays(ordinary)
            load_checkpoint!(fresh,destination,path)
            @test isequal(snapshot(fresh,destination),snapshot(st,ps))
            @test all(pair -> pair[1]===pair[2],zip(ordinaryids,H._yl_particle_arrays(ordinary)))
            @test all(candidate -> all(a -> a isa Vector, H._yl_particle_arrays(candidate)),
                fresh._stage[].candidates)
            for bad in (merge(saved,(schema_version=1,)),merge(saved,(c=nextfloat(saved.c),)),
                merge(saved,(force_policy=:rounded_gather,)),
                merge(saved,(force_policy_version=2,)),
                merge(saved,(force_policy_version=true,)),
                merge(saved,(force_policy_version=UInt(1),)),
                (; (key => value for (key,value) in pairs(saved) if key !== :force_policy)...),
                merge(saved,(last_physical_energy_change=T(Inf),)),
                merge(saved,(right_plate_charge=saved.right_plate_charge+1,)),
                merge(saved,(c_type=T===Float32 ? Float64 : Float32,)),
                merge(saved,(current_policy_version=2,)),
                merge(saved,(shape_width_bits=saved.shape_width_bits+one(saved.shape_width_bits),)),
                merge(saved,(last_dt=maxwell_cfl(restored),)),
                merge(saved,(particle_layout=collection ? :scalar : :collection,)),
                merge(saved,(fields=(E=saved.fields.E,B=(saved.fields.B[1],saved.fields.B[2],saved.fields.B[2])),)))
                before = snapshot(restored,other)
                @test_throws ArgumentError H._restore_checkpoint_state!(restored,other,bad)
                @test isequal(snapshot(restored,other),before)
            end
            step!(st,ps,T(0.002)); step!(restored,other,T(0.002))
            @test isequal(snapshot(restored,other),snapshot(st,ps))
        end
        P = typeof(ps)
        for (f,args) in ((init!,Tuple{typeof(st),P}),(step!,Tuple{typeof(st),P,T}),
            (H._checkpoint_state,Tuple{typeof(st),P}))
            report = JET.report_call(f,args;sourceinfo=:full,stacktrace_types_limit=0)
            if !isempty(JET.get_reports(report))
                show(IOContext(stdout,:limit=>false),MIME"text/plain"(),report)
            end
            @test isempty(JET.get_reports(report))
        end
    end
end
end
