"""
    CUDAHybridStepper <: AbstractIntegrator

Backend-free dispatch boundary for the optional resident CUDA D1/3V hybrid.
Construct it with `cuda_hybrid_stepper` after importing CUDA. The extension owns
the concrete storage, kernels and plans; ordinary `HybridStepper` is unchanged.
"""
abstract type CUDAHybridStepper <: AbstractIntegrator end

cuda_hybrid_stepper(args...; kwargs...) = _extension_missing(:cuda, :cuda_hybrid_stepper)
cuda_hybrid_host_snapshot(args...; kwargs...) =
    _extension_missing(:cuda, :cuda_hybrid_host_snapshot)
cuda_hybrid_host_diagnostics(args...; kwargs...) =
    _extension_missing(:cuda, :cuda_hybrid_host_diagnostics)
_validate_checkpoint_state(s, st::CUDAHybridStepper) =
    _extension_missing(:cuda, :_validate_checkpoint_state)
