# Compile-time scalar arithmetic policy. The default keeps existing promotion
# and binary operators; device callers opt in without replacing global methods.
abstract type _ScalarOps end
struct _DefaultScalarOps <: _ScalarOps end

@inline _scalar_add(::_DefaultScalarOps, a, b) = a + b
@inline _scalar_sub(::_DefaultScalarOps, a, b) = a - b
@inline _scalar_mul(::_DefaultScalarOps, a, b) = a * b
@inline _scalar_div(::_DefaultScalarOps, a, b) = a / b
