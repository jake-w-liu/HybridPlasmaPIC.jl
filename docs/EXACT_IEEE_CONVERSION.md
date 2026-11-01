# Exact IEEE conversion

Cold numerical fallbacks share the private owners in
`src/utils/exact_arithmetic.jl`. Float16/32/64 input extraction reads the stored
significand and exponent with integer arithmetic. Rounded output conversion
assembles the IEEE sign, exponent and fraction bits directly, including signed
zero, subnormal boundaries, binade carries and infinity endpoints.

Integer quotient/remainder and exact-comparison owners choose rounding before
packing. Packing performs no floating conversion or rescaling, so it preserves
that decision under hardware rounding and flush-to-zero settings. The square-root
owner retains its explicitly requested rounding mode. BigFloat keeps its
precision, rounding and unbounded-significand path.

`test/test_exact_ieee_packing.jl` contains the focused boundary and context
regressions; the complete package and compiler checks run in GitHub Actions. This
contract covers the shared conversion owners, not every compound caller or GPU
operation. Caller-specific range, accuracy and hardware contracts still apply.
