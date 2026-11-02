# Included inside HybridPlasmaPICCUDAExt, after the ordinary CUDA hooks.
const _CH = HybridPlasmaPIC
const _CHF = CUDA.CuArray{Float64,1,CUDA.DeviceMemory}
const _CHI = CUDA.CuArray{UInt64,1,CUDA.DeviceMemory}
const _CHT = CUDA.CuArray{UInt32,1,CUDA.DeviceMemory}
const _CHC = CUDA.CuArray{ComplexF64,1,CUDA.DeviceMemory}
const _CHStatus = CUDA.CuArray{Int64,1,CUDA.DeviceMemory}
const _CHParticle = _CH.ParticleSet{1,Float64,_CHF,_CHF,_CHF,_CHI,_CHT}
const _CHModel = _CH.HybridModel{_CH.IsothermalElectrons{Float64}}
const _CH_POLICY = 1
const _CH_INDEX_LIMIT = Int(typemax(Int32)) - 2

struct _CUDAHybridFields
    n::_CHF
    ui::NTuple{3,_CHF}
    B::NTuple{3,_CHF}
    E::NTuple{3,_CHF}
    J::NTuple{3,_CHF}
    pe::_CHF
    gradp::Tuple{_CHF}
    ninv::_CHF
end

_ch_field_arrays(f::_CUDAHybridFields) =
    (f.n, f.ui..., f.B..., f.E..., f.J..., f.pe, f.gradp..., f.ninv)

struct _CUDAHybridMotion
    x::Tuple{_CHF}
    v::NTuple{3,_CHF}
    xmid::Tuple{_CHF}
    vpred::NTuple{3,_CHF}
    Ep::NTuple{3,_CHF}
    Bp::NTuple{3,_CHF}
    scales::NTuple{4,_CHF}
end

_ch_motion_arrays(w::_CUDAHybridMotion) =
    (w.x..., w.v..., w.xmid..., w.vpred..., w.Ep..., w.Bp..., w.scales...)

struct _CUDAHybridGrid{M}
    n::Tuple{Int}
    L::Tuple{Float64}
    dx::Tuple{Float64}
    metadata::M
    ik::_CHC
    gain::Float64
end

struct _CUDAHybridFourier{P,Q}
    input::_CHC
    accumulator::_CHC
    forward::P
    inverse::Q
end

struct _CUDAHybridWork
    raw::NTuple{4,_CHF}
    frozen_n::_CHF
    frozen_ui::NTuple{3,_CHF}
    pe::_CHF
    gradp::Tuple{_CHF}
    ninv::_CHF
    trial::NTuple{3,_CHF}
    E::NTuple{3,_CHF}
    J::NTuple{3,_CHF}
    rates::NTuple{4,NTuple{3,_CHF}}
end

_ch_work_arrays(w::_CUDAHybridWork) = (
    w.raw...,
    w.frozen_n,
    w.frozen_ui...,
    w.pe,
    w.gradp...,
    w.ninv,
    w.trial...,
    w.E...,
    w.J...,
    w.rates[1]...,
    w.rates[2]...,
    w.rates[3]...,
    w.rates[4]...,
)

mutable struct _ResidentCUDAHybrid{G,F,C} <: _CH.CUDAHybridStepper
    g::G
    model::_CHModel
    shape::_CH.CIC
    dt::Float64
    NB::Int
    configuration::C
    species::Vector{_CHParticle}
    scalars::Vector{Tuple{Float64,Float64}}
    factors::Vector{_CH._BorisSpeciesFactors{Float64}}
    counts::Vector{Int}
    binding::Vector{Tuple{_CHF,_CHF,_CHF,_CHF,_CHF,_CHI,_CHT}}
    fields::_CUDAHybridFields
    candidate::_CUDAHybridFields
    motion::Vector{_CUDAHybridMotion}
    work::_CUDAHybridWork
    fourier::F
    status::_CHStatus
    context::CUDA.CuContext
    device::CUDA.CuDevice
    stream::CUDA.CuStream
    initialized::Bool
    primed::Bool
    time::Float64
    step::Int
    last_dt::Float64
    busy::Threads.Atomic{Bool}
end

function _ch_configuration(g, model, h, NB, scalars)
    return (
        n = g.n,
        L = g.L,
        dx = g.dx,
        gain = g.gain,
        Te = model.closure.Te,
        nfloor = model.nfloor,
        eta = model.η,
        etaH = model.ηH,
        de2 = model.de2,
        dt = h,
        NB = NB,
        shape = :CIC,
        policy = _CH_POLICY,
        scalars = Tuple(scalars),
        spectral = :release018_all_three_centered,
        deposit = :atomic_float64,
        arithmetic = :separate_rne_native_reject,
        boundary = :periodic_constant_Bx,
    )
end

_ch_particle_leaves(ps::_CHParticle) = (ps.x..., ps.v..., ps.weight, ps.id, ps.tag)
_ch_binding(ps::_CHParticle) = (ps.x[1], ps.v..., ps.weight, ps.id, ps.tag)

_ch_float(n) = CUDA.CuArray{Float64,1,CUDA.DeviceMemory}(undef, n)
_ch_vector3(n) = ntuple(_ -> _ch_float(n), 3)

function _ch_fields(n)
    return _CUDAHybridFields(
        _ch_float(n),
        _ch_vector3(n),
        _ch_vector3(n),
        _ch_vector3(n),
        _ch_vector3(n),
        _ch_float(n),
        (_ch_float(n),),
        _ch_float(n),
    )
end

function _ch_motion(n)
    return _CUDAHybridMotion(
        (_ch_float(n),),
        _ch_vector3(n),
        (_ch_float(n),),
        _ch_vector3(n),
        _ch_vector3(n),
        _ch_vector3(n),
        ntuple(_ -> _ch_float(n), 4),
    )
end

function _ch_work(n)
    return _CUDAHybridWork(
        ntuple(_ -> _ch_float(n), 4),
        _ch_float(n),
        _ch_vector3(n),
        _ch_float(n),
        (_ch_float(n),),
        _ch_float(n),
        _ch_vector3(n),
        _ch_vector3(n),
        _ch_vector3(n),
        ntuple(_ -> _ch_vector3(n), 4),
    )
end
