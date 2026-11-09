# Physical-gradient fluxes are observed separately from weak nodal reactions.
function _ace_flux(state::_AdaptiveCompositePhysical{D}) where {D}
    mesh,raw = state.mesh,_acq_raw_potential(state.mesh,state.phi)
    exterior = _CompositeQ[]
    jumps = _CompositeQ[]
    boundary_faces, interior_faces = Int[], Int[]
    for (index,face) in enumerate(mesh.faces)
        x = ntuple(d -> (face.lower[d]+face.upper[d])/2,D)
        if iszero(face.left)
            E = _acq_field(mesh,raw,face.right,x)
            push!(exterior,-face.area*E[face.axis])
            push!(boundary_faces,index)
        elseif iszero(face.right)
            E = _acq_field(mesh,raw,face.left,x)
            push!(exterior,face.area*E[face.axis])
            push!(boundary_faces,index)
        else
            left,right = _acq_field(mesh,raw,face.left,x),_acq_field(mesh,raw,face.right,x)
            push!(jumps,face.area*(left[face.axis]-right[face.axis]))
            push!(interior_faces,index)
        end
    end
    return (boundary_faces = boundary_faces,interior_faces = interior_faces,outward = exterior,normal_jumps = jumps,total_outward = sum(exterior; init = zero(_CompositeQ)))
end
"Actual piecewise Q1 gradient at a closed-box point; outer faces use interior one-sided traces."
function composite_electric_field(st::AdaptiveCompositeElectrostaticPIC{D,T},x::NTuple{D,T}) where {D,T}
    _ace_acquire(st)
    try
        _ace_ready(st,st._owners)
        all(isfinite,x) || throw(ArgumentError("composite field query requires finite coordinates"))
        exact = map(_exact_rational,x)
        k = _camr_locate(st.mesh,exact; trace = true)
        return _acq_field(st.mesh,_acq_raw_potential(st.mesh,st._state.phi),k,exact)
    finally
        _ace_release(st)
    end
end
"Detached patch volume means, including covered cells; particle forces use the full gradient."
function patch_electric_field(st::AdaptiveCompositeElectrostaticPIC{D,T},patch::Integer) where {D,T}
    _ace_acquire(st)
    try
        _ace_ready(st,st._owners)
        1 <= patch <= length(st.mesh.patches) || throw(ArgumentError("composite field patch index is invalid"))
        p = st.mesh.patches[Int(patch)]
        n = Int(patch) == 1 ? st.mesh.nbase : ntuple(_ -> 2,D)
        result = ntuple(_ -> Array{T,D}(undef,n),3)
        denominators = _aca_denominators(st.mesh.nbase,p.level)
        raw = _acq_raw_potential(st.mesh,st._state.phi)
        for (slot,a) in enumerate(p.cells)
            lo,hi = _camr_bounds(a,st.mesh.L,denominators)
            mean = _acq_box_field(st.mesh,raw,lo,hi)
            for c = 1:3
                result[c][slot] = _from_exact(T,mean[c])
                isfinite(result[c][slot]) || throw(ArgumentError("composite patch mean is nonfinite"))
            end
        end
        return result
    finally
        _ace_release(st)
    end
end
"Detached adaptive composite charge, physical-gradient/weak-Gauss, energy and full history diagnostics."
function composite_electrostatic_diagnostics(st::AdaptiveCompositeElectrostaticPIC,species::AbstractVector)
    _ace_acquire(st)
    try
        source = _ace_ready(st,species)
        state = st._state
        reaction,residual = _acq_residual(state.operator,state.phi,source.total)
        inventory = _ace_inventory(state,source)
        flux = _ace_flux(state)
        history = state.step_report === nothing ? st._history : (epoch_start_time = st.epoch_start_time,epoch_start_step = st.epoch_start_step,global_time = st.time,global_step = st.step,last_dt = st.last_dt,state = _ace_physical_payload(state))
        K,U = _kinetic_energy_rational(state.species),_acq_energy(state.operator,state.phi)
        return deepcopy(merge(_ace_phase(st),(mesh = st.mesh,phi = st.phi,E_mean = _acq_means(st.mesh,state.phi),nodal_source = source.total,inside_nodal_charge = source.inside,tail_nodal_charge = source.tail,background_nodal_charge = source.background,inventory = inventory,stratum_axes = Tuple(_acc_strata(Val(length(st.mesh.L)))),weak_gauss_residual = residual,free_conductor_charge = reaction[st.mesh.boundary],net_wall_charge = reaction[st.mesh.boundary]+source.tail[st.mesh.boundary],boundary_faces = flux.boundary_faces,interior_faces = flux.interior_faces,outward_boundary_flux = flux.outward,interior_normal_flux_jumps = flux.normal_jumps,physical_flux_reaction_defect = flux.total_outward+inventory.net_wall,kinetic_energy = K,field_energy = U,total_energy = K+U,accounting_defect = K+U-st.initial_energy-st.physical_defect-st.mesh_work,solve_report = state.solve_report,last_step_report = state.step_report,last_physical_step = history,epoch_steps = st.step-st.epoch_start_step)))
    finally
        _ace_release(st)
    end
end
