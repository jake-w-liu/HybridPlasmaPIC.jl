# Range profiles at every size and boundary mix: multipliers outside the
# normal range (Float64 power-of-two scaling, Float32 wide arithmetic) and
# overflowing native entries (exact completion) must meet the same normwise
# stage contract as ordinary stages, checked with the Rational oracle above.
function scaled_stage_fixture(::Type{T}, D, boundary, profile, n) where {T}
    wide = T === Float32 ? T(1e20) : T(1e200)
    small = inv(wide)
    speed = profile === :c_large ? wide : profile === :c_small ? small : one(T)
    L = profile === :derivative_large ? small : profile === :derivative_small ? wide : one(T)
    bc = ntuple(d -> boundary === :mixed ? (isodd(d) ? :pec : :periodic) : boundary, D)
    m = MaxwellVacuum(YeeGrid(ntuple(_ -> n, D), ntuple(_ -> L, D); boundaries = bc); c = speed)
    dt = profile === :half_subnormal ? nextfloat(zero(T)) : T(0.1) * maxwell_cfl(m)
    amplitude =
        profile === :c_large ? small :
        profile === :c_small ? wide :
        profile === :derivative_large ? wide :
        profile === :derivative_small ? small : floatmax(T) / 16
    for c = 1:3, I in CartesianIndices(m.fields.B[c])
        m.fields.B[c][I] = T(mod(sum(Tuple(I)) + c, 3) - 1) * amplitude
    end
    if profile in (:derivative_large, :derivative_small, :half_subnormal)
        for c = 1:3, I in CartesianIndices(m.fields.E[c])
            m.fields.E[c][I] = T(mod(sum(Tuple(I)) + c, 3) - 1) * amplitude
        end
    end
    apply_yee_pec!(m.fields, m.g)
    m, dt
end

@testset "Range-profile stages meet the normwise contract" begin
    for T in (Float32, Float64),
        D = 1:3,
        boundary in (D == 1 ? (:periodic, :pec) : (:periodic, :pec, :mixed)),
        profile in (:c_large, :c_small, :derivative_large, :derivative_small, :half_subnormal),
        n in (2, 4, 8)

        m, dt = scaled_stage_fixture(T, D, boundary, profile, n)
        normwise_step_check(m, dt)
        @test m.time[] == dt && m.step[] == 1
    end
end
