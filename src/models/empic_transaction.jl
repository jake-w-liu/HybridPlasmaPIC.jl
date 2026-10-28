# Undo storage for public relativistic EMPIC steps. The live solver keeps its
# existing arrays and FFT plans; only their enumerated values are copied.
mutable struct _EMPICSpeciesSnapshot{D,T}
    x::NTuple{D,Vector{T}}
    v::NTuple{3,Vector{T}}
    weight::Vector{T}
    id::Vector{UInt64}
    tag::Vector{UInt32}
    q::T
    m::T
end

mutable struct _EMPICTransaction{D,T}
    electrons::_EMPICSpeciesSnapshot{D,T}
    ions::Union{Nothing,_EMPICSpeciesSnapshot{D,T}}
    fields::Vector{Array{T,D}}
    electron_gathers::Vector{Vector{T}}
    ion_gathers::Vector{Vector{T}}
    fft::NTuple{4,Array{Complex{T},D}}
    clock::Tuple{T,Int,T,T}
end

function _EMPICSpeciesSnapshot(ps::ParticleSet{D,T}) where {D,T}
    return _EMPICSpeciesSnapshot{D,T}(
        ntuple(d -> Vector{T}(undef, length(ps.x[d])), Val(D)),
        ntuple(d -> Vector{T}(undef, length(ps.v[d])), Val(3)),
        Vector{T}(undef, length(ps.weight)),
        Vector{UInt64}(undef, length(ps.id)),
        Vector{UInt32}(undef, length(ps.tag)),
        ps.q,
        ps.m,
    )
end

@inline function _empic_snapshot_vector!(saved::Vector, live::AbstractVector)
    length(saved) == length(live) || resize!(saved, length(live))
    copyto!(saved, live)
    return nothing
end

function _empic_save_species!(saved::_EMPICSpeciesSnapshot{D}, live::ParticleSet{D}) where {D}
    for d = 1:D
        _empic_snapshot_vector!(saved.x[d], live.x[d])
    end
    for d = 1:3
        _empic_snapshot_vector!(saved.v[d], live.v[d])
    end
    _empic_snapshot_vector!(saved.weight, live.weight)
    _empic_snapshot_vector!(saved.id, live.id)
    _empic_snapshot_vector!(saved.tag, live.tag)
    saved.q, saved.m = live.q, live.m
    return nothing
end

function _empic_restore_species!(live::ParticleSet{D}, saved::_EMPICSpeciesSnapshot{D}) where {D}
    for d = 1:D
        copyto!(live.x[d], saved.x[d])
    end
    for d = 1:3
        copyto!(live.v[d], saved.v[d])
    end
    copyto!(live.weight, saved.weight)
    copyto!(live.id, saved.id)
    copyto!(live.tag, saved.tag)
    live.q, live.m = saved.q, saved.m
    return nothing
end

function _new_empic_transaction(es, electrons::ParticleSet{D,T}, ions) where {D,T}
    return _EMPICTransaction{D,T}(
        _EMPICSpeciesSnapshot(electrons),
        ions === nothing ? nothing : _EMPICSpeciesSnapshot(ions),
        Array{T,D}[similar(a) for a in _empic_transaction_fields(es)],
        Vector{T}[similar(a) for a in _empic_electron_gathers(es)],
        Vector{T}[similar(a) for a in _empic_ion_gathers(es)],
        map(similar, (es.g.cbuf, es.g.tbuf, es.g.abuf, es.g.sbuf)),
        (es.time[], es.step[], es.last_dt[], es.last_sub_dt[]),
    )
end

function _save_empic_transaction!(es, electrons, ions)
    saved = es.transaction
    if saved === nothing
        saved = _new_empic_transaction(es, electrons, ions)
        es.transaction = saved
    end
    _empic_save_species!(saved.electrons, electrons)
    if ions !== nothing
        saved_ions = saved.ions
        if saved_ions === nothing
            saved_ions = _EMPICSpeciesSnapshot(ions)
            saved.ions = saved_ions
        end
        _empic_save_species!(saved_ions, ions)
    end
    fields = _empic_transaction_fields(es)
    for i in eachindex(saved.fields)
        copyto!(saved.fields[i], fields[i])
    end
    for (backups, live) in (
        (saved.electron_gathers, _empic_electron_gathers(es)),
        (saved.ion_gathers, _empic_ion_gathers(es)),
    )
        for i in eachindex(backups)
            _empic_snapshot_vector!(backups[i], live[i])
        end
    end
    buffers = (es.g.cbuf, es.g.tbuf, es.g.abuf, es.g.sbuf)
    for i = 1:4
        copyto!(saved.fft[i], buffers[i])
    end
    saved.clock = (es.time[], es.step[], es.last_dt[], es.last_sub_dt[])
    return saved
end

function _restore_empic_transaction!(es, electrons, ions, saved)
    # Restore original gather lengths as well as values. These buffers may have
    # been resized after the snapshot to accommodate a changed population.
    for (live, backups) in (
        (_empic_electron_gathers(es), saved.electron_gathers),
        (_empic_ion_gathers(es), saved.ion_gathers),
    )
        for i in eachindex(backups)
            _empic_snapshot_vector!(live[i], backups[i])
        end
    end
    fields = _empic_transaction_fields(es)
    for i in eachindex(saved.fields)
        copyto!(fields[i], saved.fields[i])
    end
    buffers = (es.g.cbuf, es.g.tbuf, es.g.abuf, es.g.sbuf)
    for i = 1:4
        copyto!(buffers[i], saved.fft[i])
    end
    _empic_restore_species!(electrons, saved.electrons)
    if ions !== nothing
        saved_ions = saved.ions
        saved_ions === nothing && error("missing EMPIC ion transaction snapshot")
        _empic_restore_species!(ions, saved_ions)
    end
    es.time[], es.step[], es.last_dt[], es.last_sub_dt[] = saved.clock
    return nothing
end

function _empic_step_transaction!(es, electrons, ions, dt, sub_dt, ns, next_time, next_step)
    required_ions = es.mobile ? ions : nothing
    saved = _save_empic_transaction!(es, electrons, required_ions)
    completed = false
    try
        _advance_empic_step!(es, electrons, ions, dt, sub_dt, ns, next_time, next_step)
        completed = true
        return es
    finally
        completed || _restore_empic_transaction!(es, electrons, required_ions, saved)
    end
end
