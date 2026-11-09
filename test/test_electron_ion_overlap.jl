module ElectronIonOverlapTests
using HybridPlasmaPIC, Test, Random, Serialization
const H = HybridPlasmaPIC
const Q = Rational{BigInt}

function fixture(::Type{T}, shape; D=1, dims=ntuple(d -> d == 1 ? 4 : 2,D),
    ns=2, cold=false, sparse=false, enabled=true, logs=one(T),
    pressure=T(4), bulk=zero(T), epsilon=T(2)^(-28), R=Xoshiro,
    closure=EvolvingElectrons(T(5/3))) where {T}
    g = FourierGrid(dims,Tuple(T.(dims)))
    species = [ParticleSet{D,T}(3prod(dims);q=T(a),m=T(a+1)) for a in 1:ns]
    for (a,p) in enumerate(species), I in CartesianIndices(dims), j in 1:3
        frame = LinearIndices(dims)[I]
        i = 3(frame-1)+j
        for d in 1:D
            p.x[d][i] = T(I[d]-1)+T(3/8)
        end
        p.weight[i] = j == 1 ? one(T) : j == 2 && !(sparse && frame == prod(dims)) ? T(3) : zero(T)
        amplitude = T(a+mod(frame,3))/T(32)
        for c in 1:3
            p.v[c][i] = j == 3 ? zero(T) : bulk + (cold ? zero(T) :
                (j == 1 ? T(3) : -one(T))*amplitude*T((1,-1/2,1/4)[c]))
        end
    end
    inner = HybridStepper(g,HybridModel(closure),shape,sum(nparticles,species))
    fill!(inner.fields.pe,pressure)
    request = enabled ? ElectronIonThermalExchange(;electron_mass=T(1/16),coulomb_log=logs) : nothing
    # This exact interval keeps the dyadic off-node positions fixed in the
    # ordinary fixtures. Small epsilon keeps the source transfer resolvable.
    st = CollisionalHybridStepper(inner,species;dt=T(2)^(-60),rng=R(913),
        coulomb_log=zero(T),epsilon0=epsilon,weight_to_number=one(T),electron_ion=request)
    init!(st,species)
    return st,species
end

function bytes(value)
    io = IOBuffer()
    serialize(io,value)
    return take!(io)
end
function snapshot(st,species)
    s = H._checkpoint_state(st,species)
    return deepcopy(merge(s,(rng=bytes(s.rng),)))
end
function frameof(p,g,i)
    I = ntuple(d -> mod(floor(Int,Q(p.x[d][i])/Q(g.dx[d])),g.n[d])+1,length(g.n))
    return LinearIndices(g.n)[I...]
end
function moments(p,g,frame)
    ids = [i for i in eachindex(p.weight) if p.weight[i]>0 && frameof(p,g,i)==frame]
    W,m = sum((Q(p.weight[i]) for i in ids);init=zero(Q)),Q(p.m)
    P = ntuple(c -> sum((m*Q(p.weight[i])*Q(p.v[c][i]) for i in ids);init=zero(Q)),3)
    mean = iszero(W) ? ntuple(_ -> zero(Q),3) : ntuple(c -> P[c]/(m*W),3)
    K = sum((m*Q(p.weight[i])*sum(Q(p.v[c][i])^2 for c in 1:3)/2 for i in ids);init=zero(Q))
    U = sum((m*Q(p.weight[i])*sum((Q(p.v[c][i])-mean[c])^2 for c in 1:3)/2 for i in ids);init=zero(Q))
    return (weight=W,mean=mean,kinetic=K,momentum=P,thermal=U,capacity=3W/2,count=length(ids))
end

# Independent cardinal polynomials at dyadic coordinates. Each axis weight is
# rounded separately, then the stored factors are multiplied as exact Q values.
# No production stencil, incidence visitor, bin or rate helper is called.
function incidence(species,g,shape)
    out = Dict{NTuple{3,Int},Q}()
    for (a,p) in enumerate(species), i in eachindex(p.weight)
        iszero(p.weight[i]) && continue
        stencils = ntuple(length(g.n)) do d
            T = eltype(p.weight)
            s = Q(p.x[d][i])/Q(g.dx[d])
            c = shape isa CIC ? floor(Int,s) : round(Int,s)
            f = s-c
            weights = shape isa CIC ? (1-f,f) : ((1//2)*(1//2-f)^2,3//4-f^2,(1//2)*(1//2+f)^2)
            return (base=shape isa CIC ? c : c-1,weights=map(x -> Q(T(x)),weights))
        end
        for O in CartesianIndices(ntuple(d -> length(stencils[d].weights),length(g.n)))
            I = ntuple(d -> mod(stencils[d].base+O[d]-1,g.n[d])+1,length(g.n))
            weight = Q(p.weight[i])*prod(stencils[d].weights[O[d]] for d in 1:length(g.n))
            iszero(weight) && continue
            key = (a,frameof(p,g,i),LinearIndices(g.n)[I...])
            out[key] = get(out,key,zero(Q))+weight
        end
    end
    return [(species=a,frame=f,node=j,weight=out[(a,f,j)]) for (a,f,j) in sort!(collect(keys(out)))]
end

function physical_beta(st,edge,initial,density,pressure)
    return setprecision(BigFloat,256) do
        p = st.species[edge.species]
        V = prod(BigFloat,st.grid.dx)
        Ce = BigFloat(density[edge.node])*V/(BigFloat(st.inner.model.closure.γ)-1)
        Te = BigFloat(pressure[edge.node])/BigFloat(density[edge.node])
        Ti = BigFloat(initial.thermal)/BigFloat(initial.capacity)
        me,mi = BigFloat(st.electron_ion.electron_mass),BigFloat(p.m)
        rate = 8sqrt(2BigFloat(pi))*BigFloat(st.weight_to_number)*BigFloat(edge.weight)/V*
            BigFloat(p.q)^2*BigFloat(st.electron_ion.logarithm[edge.species]) /
            (3*(4BigFloat(pi)*BigFloat(st.epsilon0))^2*me*mi*(Te/me+Ti/mi)^(BigFloat(3)/2))
        lambda = BigFloat(st.dt)*Ce*rate/BigFloat(initial.capacity)
        return lambda/(1+lambda)
    end
end

function check_report(st,species,reference,refspecies)
    r = electron_ion_exchange_diagnostics(st,species)
    T,g = eltype(st.inner.fields.pe),st.grid
    nc,V,gamma = prod(g.n),prod(Q,g.dx),Q(st.inner.model.closure.γ)-1
    @test r.policy == 2 && r.frame_policy === :periodic_floor_cell && r.edge_order === :species_frame_node
    @test r.shape === (st.inner.shape isa CIC ? :cic : :tsc)
    @test all(a -> species[a].x == refspecies[a].x,eachindex(species))
    @test r.initial_pressure == reference.inner.fields.pe
    @test r.density == st.inner.fields.n
    @test r.incidence == incidence(refspecies,g,st.inner.shape)
    initial = [moments(p,g,f) for p in refspecies for f in 1:nc]
    @test r.initial_frames == [(species=a,frame=f,statistics=initial[(a-1)*nc+f]) for a in eachindex(species) for f in 1:nc]
    partition = Q[-s.weight for s in initial]
    charge = map(x -> Q(x)*V,r.density)
    for edge in r.incidence
        partition[(edge.species-1)*nc+edge.frame] += edge.weight
        charge[edge.node] -= Q(species[edge.species].q)*edge.weight
    end
    @test r.frame_partition == partition && r.charge_defect == charge
    # One exact incidence serves both directions of this adjoint identity.
    electron = Q[Q(j+1)/7 for j in 1:nc]
    ion = Q[Q(f+2)/11 for f in eachindex(initial)]
    at_nodes,at_frames = zeros(Q,nc),zeros(Q,length(initial))
    for e in r.incidence
        f = (e.species-1)*nc+e.frame
        at_nodes[e.node] += e.weight*ion[f]
        at_frames[f] += e.weight*electron[e.node]
    end
    @test sum(electron.*at_nodes) == sum(ion.*at_frames)
    @test [(e.species,e.frame,e.node,e.weight) for e in r.edges] ==
        [(e.species,e.frame,e.node,e.weight) for e in r.incidence if st.electron_ion.logarithm[e.species]>0]
    current,pressure = copy(initial),copy(r.initial_pressure)
    total_defect = total_bound = zero(Q)
    momentum_defect = momentum_bound = ntuple(_ -> zero(Q),3)
    for e in r.edges
        f,j = (e.species-1)*nc+e.frame,e.node
        @test e.before == current[f] && isequal(e.pressure_before,pressure[j])
        b = physical_beta(st,e,initial[f],r.density,r.initial_pressure)
        @test e.beta[1] <= Q(b) <= e.beta[3]
        @test e.beta[1] <= e.beta[2] <= e.beta[3] < 1
        Ce,Te,Ci,Ti = Q(r.density[j])*V/gamma,Q(e.pressure_before)/Q(r.density[j]),e.before.capacity,e.before.thermal/e.before.capacity
        conductance_dt = Ci*e.beta[2]/(1-e.beta[2])
        heat = iszero(conductance_dt) ? zero(Q) : (Te-Ti)/(inv(conductance_dt)+inv(Ce)+inv(Ci))
        @test e.target == e.before.thermal+heat
        @test e.electron_target == Te-heat/Ce
        deltaK = e.after.kinetic-e.before.kinetic
        defect = deltaK+V*(Q(e.pressure_after)-Q(e.pressure_before))/gamma
        @test defect == e.energy_defect && abs(defect) <= e.energy_bound
        @test all(c -> abs(e.after.momentum[c]-e.before.momentum[c]) <= e.materialized.momentum_bound[c],1:3)
        @test abs(e.after.kinetic-e.materialized.frame_energy) <= e.materialized.kinetic_bound
        @test e.entropy in (:resolved,:cold_positive,:unchanged)
        @test e.entropy_bounds[1] >= 0
        # Independently bound historical tight certificates by IEEE spacing
        # and weighted Cauchy, using directed high-precision upper roots.
        mp_cap,k_cap = setprecision(BigFloat,512) do
            setrounding(BigFloat,RoundUp) do
                m,W,K = BigFloat(species[e.species].m),BigFloat(e.after.weight),BigFloat(e.after.kinetic)
                u,eta = BigFloat(2)^(-precision(T)),BigFloat(nextfloat(zero(T)))
                S = 2K/m
                (Q(m*(u*sqrt(W*S)+eta*W)),Q(m*((u+3u^2/2)*S+(1+3u)*eta*sqrt(3W*S)+9eta^2*W/2)))
            end
        end
        @test all(x -> x <= mp_cap,e.materialized.momentum_bound)
        @test e.materialized.kinetic_bound <= k_cap
        total_defect += defect
        total_bound += e.energy_bound
        momentum_defect = ntuple(c -> momentum_defect[c]+e.after.momentum[c]-e.before.momentum[c],3)
        momentum_bound = ntuple(c -> momentum_bound[c]+e.materialized.momentum_bound[c],3)
        current[f],pressure[j] = e.after,e.pressure_after
    end
    @test current == [moments(p,g,f) for p in species for f in 1:nc]
    @test pressure == st.inner.fields.pe
    energy_change = sum(current[f].kinetic-initial[f].kinetic for f in eachindex(initial))+
        V*sum(Q(pressure[j])-Q(r.initial_pressure[j]) for j in eachindex(pressure))/gamma
    @test energy_change == r.energy_defect == total_defect
    @test abs(total_defect) <= r.energy_bound == total_bound
    @test r.momentum_defect == momentum_defect && r.momentum_bound == momentum_bound
    @test all(c -> abs(momentum_defect[c]) <= momentum_bound[c],1:3)
    return r
end

@testset "CIC/TSC stored overlap, physical rates, causal energy and actual feedback" begin
    for T in (Float32,Float64), shape in (CIC(),TSC()), D in 1:3
        st,ps = fixture(T,shape;D)
        ref,rps = fixture(T,shape;D,enabled=false)
        step!(ref,rps,ref.dt)
        step!(st,ps,st.dt)
        report = check_report(st,ps,ref,rps)
        @test any(abs.(st.inner.fields.pe-ref.inner.fields.pe).>0)
        @test any(maximum(abs,e)>0 for e in st.inner.fields.E)
        detached = electron_ion_exchange_diagnostics(st,ps)
        empty!(detached.edges)
        @test !isempty(electron_ion_exchange_diagnostics(st,ps).edges)
        @test all(a -> all(c -> isequal(ps[a].v[c][3:3:end],rps[a].v[c][3:3:end]),1:3),eachindex(ps))
    end
    # TSC's 27 images collapse onto two actual nodes on this periodic mesh.
    st,ps = fixture(Float64,TSC();D=3,dims=(2,1,1))
    ref,rps = fixture(Float64,TSC();D=3,dims=(2,1,1),enabled=false)
    step!(ref,rps,ref.dt); step!(st,ps,st.dt)
    r = check_report(st,ps,ref,rps)
    @test length(r.incidence) == 8
    # F32 rounds 1-2^-26 to one. Keep its physical frame capacity distinct
    # from the resulting stored-shape sum, rather than renormalizing H.
    st,ps = fixture(Float32,CIC();dims=(2,),ns=1)
    ref,rps = fixture(Float32,CIC();dims=(2,),ns=1,enabled=false)
    for p in (ps[1],rps[1])
        p.x[1][1:3] .= Float32(2)^(-26)
    end
    init!(st,ps); init!(ref,rps)
    step!(ref,rps,ref.dt); step!(st,ps,st.dt)
    r = check_report(st,ps,ref,rps)
    @test r.frame_partition[1] == 4Q(Float32(2)^(-26)) > 0
    st,ps = fixture(Float64,CIC();closure=TransportElectrons(5/3))
    ref,rps = fixture(Float64,CIC();closure=TransportElectrons(5/3),enabled=false)
    step!(ref,rps,ref.dt); step!(st,ps,st.dt)
    check_report(st,ps,ref,rps)
    gradient = similar(st.inner.fields.pe)
    deriv!(gradient,st.inner.fields.pe,st.grid,1)
    force = -gradient./st.inner.fields.n
    @test isapprox(st.inner.pressure_work.Eforce[1],force;rtol=128eps(),atol=128eps()*maximum(abs,force))
end

@testset "Overlap cold draws, global sparse preflight, late undo and zero reduction" begin
    for T in (Float32,Float64), shape in (CIC(),TSC()), R in (Xoshiro,MersenneTwister)
        st,ps = fixture(T,shape;ns=1,cold=true,R)
        expected = deepcopy(st.rng)
        for _ in 1:24
            randn(expected,T)
        end
        step!(st,ps,st.dt)
        @test bytes(st.rng) == bytes(expected)
        @test all(e -> e.after.thermal>0,electron_ion_exchange_diagnostics(st,ps).edges)
        sparse,sp = fixture(T,shape;ns=1,cold=true,sparse=true,R)
        entry = snapshot(sparse,sp)
        @test_throws ArgumentError step!(sparse,sp,sparse.dt)
        @test isequal(snapshot(sparse,sp),entry)
        # The small heat is below the large bulk ULP, after cold RNG draws.
        late,lp = fixture(T,shape;ns=1,cold=true,R,bulk=T(2)^50,pressure=T(1e-6),epsilon=one(T))
        entry = snapshot(late,lp)
        failure = try
            step!(late,lp,late.dt); nothing
        catch err
            err
        end
        @test failure isa ArgumentError
        @test occursin("lost positive thermal energy",sprint(showerror,failure))
        @test isequal(snapshot(late,lp),entry)
        for p in lp, v in p.v
            fill!(v,zero(T))
        end
        init!(late,lp)
        clean,cp = fixture(T,shape;ns=1,cold=true,R,pressure=T(1e-6),epsilon=one(T))
        step!(late,lp,late.dt); step!(clean,cp,clean.dt)
        @test isequal(snapshot(late,lp),snapshot(clean,cp))
        inactive,ip = fixture(T,shape;logs=zero(T))
        ref,rp = fixture(T,shape;enabled=false)
        entry = snapshot(inactive,ip)
        step!(inactive,ip,zero(T))
        @test isequal(snapshot(inactive,ip),entry)
        step!(inactive,ip,inactive.dt); step!(ref,rp,ref.dt)
        @test H._hybrid_field_state(inactive.inner.fields) == H._hybrid_field_state(ref.inner.fields)
        @test [H._checkpoint_particle_state(p) for p in ip] == [H._checkpoint_particle_state(p) for p in rp]
        @test bytes(inactive.rng) == bytes(ref.rng)
        @test isempty(electron_ion_exchange_diagnostics(inactive,ip).edges)
        equilibrium,ep = fixture(T,shape;ns=1,cold=true,pressure=zero(T),R)
        rng = bytes(equilibrium.rng)
        step!(equilibrium,ep,equilibrium.dt)
        @test bytes(equilibrium.rng) == rng
        @test all(e -> e.cold && e.entropy === :unchanged,
            electron_ion_exchange_diagnostics(equilibrium,ep).edges)
    end
end

@testset "Overlap schema3 raw/archive continuation and malformed history undo" begin
    st,ps = fixture(Float64,CIC();ns=1,cold=true)
    step!(st,ps,st.dt)
    saved = deepcopy(H._checkpoint_state(st,ps))
    @test saved.schema == 3 && saved.electron_ion.policy == 2
    mktempdir() do dir
        raw,archive = joinpath(dir,"overlap.chk"),joinpath(dir,"overlap.run")
        save_checkpoint(raw,st,ps); archive_run(archive,st,ps;rng_seed=913)
        target,tp = fixture(Float64,CIC();ns=1,cold=true)
        step!(target,tp,target.dt); step!(target,tp,target.dt)
        identities = [map(objectid,(p.x...,p.v...,p.weight,p.id,p.tag)) for p in tp]
        load_checkpoint!(target,tp,raw)
        @test isequal(snapshot(target,tp),snapshot(st,ps))
        @test [map(objectid,(p.x...,p.v...,p.weight,p.id,p.tag)) for p in tp] == identities
        step!(target,tp,target.dt); step!(st,ps,st.dt)
        @test isequal(snapshot(target,tp),snapshot(st,ps))
        load_archive!(target,tp,archive)
        baseline,bp = fixture(Float64,CIC();ns=1,cold=true)
        H._restore_checkpoint_state!(baseline,bp,saved)
        @test isequal(snapshot(target,tp),snapshot(baseline,bp))
        for damage in (:omitted,:reordered,:infinite,:cap,:alias,:schema)
            incoming = deepcopy(saved)
            edges = incoming.thermal_report.edges
            if damage === :omitted
                pop!(edges)
            elseif damage === :reordered
                reverse!(edges)
            elseif damage === :infinite
                edges[1] = merge(edges[1],(energy_defect=Q(1,0),))
            elseif damage === :cap
                m = edges[1].materialized
                edges[1] = merge(edges[1],(materialized=merge(m,(kinetic_bound=m.kinetic_bound+1,)),))
            elseif damage === :alias
                incoming = merge(incoming,(thermal_report=merge(incoming.thermal_report,(edges=target.thermal_report.edges,)),))
            else
                incoming = merge(incoming,(schema=2,))
            end
            entry = snapshot(target,tp)
            @test_throws ArgumentError H._restore_checkpoint_state!(target,tp,incoming)
            @test isequal(snapshot(target,tp),entry)
        end
    end
end
end # module
