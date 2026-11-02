# Real mixed geometry: the normal stencil folds on physical endpoint nodes;
# transverse coordinates and the shape arithmetic retain their canonical owner.
function _particle_cell_position(ps::ParticleSet{D,T}, op::MixedShockOperators{D,T},
    d::Int, p::Int) where {D,T}
    x = ps.x[d][p]
    isfinite(x) && (d == 1 ? 0 <= x <= op.L[d] : 0 <= x < op.L[d]) ||
        throw(ArgumentError("mixed CIC positions are outside their physical domain"))
    s = d == 1 ? x / op.dx[d] : _periodic_cell_position(x, op.L[d], op.dx[d])
    isfinite(s) && _exact_rational(typemin(Int)) <= _exact_rational(s) <=
        _exact_rational(typemax(Int) - 2) || throw(ArgumentError("mixed CIC cell position is not safely Int-convertible"))
    return s
end

_open_hall_volume(op, ix) = _exact_rational(op.sbp.H[ix]) *
    prod(_exact_rational, Base.tail(op.dx))

function _open_hall_moments!(st::OpenHallHybrid{D,T}, f::HybridFields{D,T}, species) where {D,T}
    sums = st.work.moment
    for a in sums
        fill!(a, zero(_OpenHallQ))
    end
    for ps in species
        q = _exact_rational(ps.q)
        _foreach_exact_particle_contribution(ps, (ps.weight,), st.op, CIC(), Val(:fold_x)) do idx, value, p
            charge = q * value
            sums[1][idx...] += charge
            for c = 1:3
                sums[c + 1][idx...] += charge * _exact_rational(ps.v[c][p])
            end
        end
    end
    for I in CartesianIndices(f.n)
        volume = _open_hall_volume(st.op, I[1])
        f.n[I] = _from_exact(T, sums[1][I] / volume)
        for c = 1:3
            f.ui[c][I] = _from_exact(T, sums[c + 1][I] / volume)
        end
    end
    floor = _open_hall_physics(st.model, T)
    all(x -> isfinite(x) && x >= floor, f.n) && all(a -> all(isfinite, a), f.ui) ||
        throw(ArgumentError("open Hall complete moments require finite materialized density/momentum and strict n>=nfloor"))
    _normalize_ion_moments!(f.n, f.ui, floor)
    all(a -> all(isfinite, a), f.ui) || throw(ArgumentError("open Hall normalized ion velocity crosses its materialized range"))
    return nothing
end

function _open_hall_complete_fields!(st::OpenHallHybrid, f::HybridFields)
    w = st.work.hall
    # Init/endpoint/force need the completed local field, not a differentiated
    # source certificate. The RHS separately retains its canonical risk reject.
    _mixed_hall_electric!(f.B, f.n, f.ui, st.model, st.op, w)
    copyto!(f.pe, w.pe)
    copyto!(f.ninv, w.ninv)
    for (out, input) in ((f.E, w.E), (f.J, w.J), (f.gradp, w.gradp), (f.lapJ, w.lapJ))
        for c in eachindex(out)
            copyto!(out[c], input[c])
        end
    end
    f.floor_count[] = w.floor_count[]
    return nothing
end

function _open_hall_gather!(st::OpenHallHybrid{D,T}, ps, fields) where {D,T}
    n = nparticles(ps)
    Ep, Bp = st.work.Ep, st.work.Bp
    for a in (Ep..., Bp...)
        length(a) >= n || throw(DimensionMismatch("open Hall particle gather capacity was not prepared"))
    end
    stamp = CartesianIndices(ntuple(_ -> 2, D))
    for p = 1:n
        stencil = ntuple(d -> _stencil1d(CIC(), _particle_cell_position(ps, st.op, d, p)), D)
        E = _gather_vector_stencil(stencil, st.op.n, stamp, T, Val(:fold_x), Val(true)) do idx
            ntuple(c -> fields.E[c][idx...], 3)
        end
        B = _gather_vector_stencil(stencil, st.op.n, stamp, T, Val(:fold_x), Val(true)) do idx
            ntuple(c -> fields.B[c][idx...], 3)
        end
        for c = 1:3
            Ep[c][p], Bp[c][p] = E[c], B[c]
        end
    end
    return nothing
end

function _open_hall_grid_charge(st, f)
    result = zero(_OpenHallQ)
    for I in CartesianIndices(f.n)
        result += _exact_rational(f.n[I]) * _open_hall_volume(st.op, I[1])
    end
    return result
end
