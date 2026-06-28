# PIC1D HDF5 Extension Validation

This case validates one-dimensional PIC models and the HDF5 IO extension.

## What Is Tested

- Electrostatic1D Poisson solve against an exact field.
- Electrostatic field-energy calculation.
- EMPIC1D charge-conservation residual and finite-energy behavior.
- Mobile-ion/subcycled EMPIC1D charge conservation.
- HDF5 dense-array write/read extension roundtrip.

## Reference

The references are an exact one-dimensional Fourier-mode Poisson field, direct
field-energy arithmetic, charge-conservation residual thresholds, and exact HDF5
array roundtrip equality.

## Outputs

The runner writes metric CSV rows in this folder's `artifacts/` directory. The
plotter writes a PlotlySupply PDF comparing error ratios to tolerances.

## Run

```sh
julia --project=validation validation/run_validation.jl --case 16_pic1d_hdf5_extension_validation
```
