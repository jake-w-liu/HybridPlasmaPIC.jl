using Test
using LinearAlgebra
using HybridPlasmaPIC

# The reference uses only the independently assembled physical generator.
function _implicit_yee_cayley(A, initial, h)
    Q = Rational{BigInt}
    identity = Matrix{Q}(I, length(initial), length(initial))
    half = Q(h) / 2
    return (identity - half * A) \ ((identity + half * A) * Q.(initial))
end

@testset "Independent compatible Yee transverse oracles" begin
    Q = Rational{BigInt}
    for T in (Float32, Float64)
        n, h = 4, T(1 / 8)
        g = YeeGrid((n,), (T(n),); boundaries = (:periodic,))
        @testset "vacuum Cayley / true Yee Nyquist $T" begin
            species = ParticleSet{1,T}[]
            st = ImplicitYeeEMPICStepper(
                g,
                species;
                shape = YeeTriangularCloud(one(T)),
                c = one(T),
                nonlinear_rtol = 32eps(T),
                nonlinear_atol = 32eps(T),
            )
            fill!(st.fields.B[1], -zero(T))
            st.fields.E[2] .= T.((1, -1, 1, -1)) ./ T(8)
            st.fields.E[3] .= T.((0, 1, 0, -1)) ./ T(16)
            bx = copy(st.fields.B[1])
            initial = vcat(st.fields.E[2:3]..., st.fields.B[2:3]...)
            init!(st, species)
            A = zeros(Q, 4n, 4n) # (Ey, Ez, By, Bz), dx = c = 1.
            for i = 1:n
                previous, next = mod1(i - 1, n), mod1(i + 1, n)
                # Ey'=-D- Bz, Ez'=D- By; By'=D+ Ez, Bz'=-D+ Ey.
                A[i, 3n+i], A[i, 3n+previous] = -1, 1
                A[n+i, 2n+i], A[n+i, 2n+previous] = 1, -1
                A[2n+i, n+next], A[2n+i, n+i] = 1, -1
                A[3n+i, next], A[3n+i, i] = -1, 1
            end
            reference = _implicit_yee_cayley(A, initial, h)
            step!(st, species, h)
            r = st.last_report
            endpoint = vcat(st.fields.E[2:3]..., st.fields.B[2:3]...)
            # Skew generator: inverse(I-hA/2) has norm <= 1 (16 rows).
            budget = 8 * (Q(r.tolerance) + Q(eps(T)))
            @test maximum(abs, Q.(endpoint) - reference) <= budget
            @test !all(iszero, st.fields.B[2]) && !all(iszero, st.fields.B[3])
            @test isequal(st.fields.B[1], bx)
            @test all(iszero, st.fields.E[1]) && all(a -> all(iszero, a), st.J)
            @test all(iszero, st.charge_n) && all(iszero, st.charge_np1)
            @test r.gauss <= Q(256eps(T)) && r.continuity <= Q(256eps(T))
            @test r.energy_change == sum(r.defects)
            @test iszero(r.defects[2]) && iszero(r.curl_adjoint_work)
        end
        @testset "magnetized unequal q/m cold-plasma DC midpoint $T" begin
            electrons = ParticleSet{1,T}(n; q = -one(T), m = one(T))
            ions = ParticleSet{1,T}(n; q = T(2), m = T(4))
            species = [electrons, ions]
            for ps in species
                ps.x[1] .= T.(0:(n-1))
            end
            fill!(ions.weight, T(1 / 2)) # Equal opposite charge densities, dx = 1.
            fill!(electrons.v[2], T(1 / 8))
            fill!(electrons.v[3], -T(1 / 16))
            fill!(ions.v[2], -T(1 / 16))
            fill!(ions.v[3], T(1 / 32))
            st = ImplicitYeeEMPICStepper(
                g,
                species;
                shape = YeeTriangularCloud(one(T)),
                c = one(T),
                nonlinear_rtol = 32eps(T),
                nonlinear_atol = 32eps(T),
            )
            fill!(st.fields.B[1], T(1 / 2))
            fill!(st.fields.E[2], T(1 / 16))
            fill!(st.fields.E[3], -T(1 / 32))
            positions, bx = [copy(ps.x[1]) for ps in species], copy(st.fields.B[1])
            init!(st, species)
            initial = Q.([1 // 8, -1 // 16, -1 // 16, 1 // 32, 1 // 16, -1 // 32])
            A = zeros(Q, 6, 6) # (vey, vez, viy, viz, Ey, Ez), Bx = 1/2.
            for (row, factor) in ((1, -one(Q)), (3, Q(1 // 2)))
                A[row, row+1], A[row, 5] = factor / 2, factor
                A[row+1, row], A[row+1, 6] = -factor / 2, factor
            end
            A[5, 1], A[5, 3] = 1, -1 # -J = ve - vi, original q*weight.
            A[6, 2], A[6, 4] = 1, -1
            reference = _implicit_yee_cayley(A, initial, h)
            step!(st, species, h)
            r = st.last_report
            # Positive energy weights (1,1,2,2,1,1) bound the inverse midpoint map.
            budget = 16 * (Q(r.tolerance) + Q(eps(T)))
            endpoint = (
                electrons.v[2],
                electrons.v[3],
                ions.v[2],
                ions.v[3],
                st.fields.E[2],
                st.fields.E[3],
            )
            @test all(i -> maximum(x -> abs(Q(x) - reference[i]), endpoint[i]) <= budget, 1:6)
            for (component, electron, ion) in ((2, 1, 3), (3, 2, 4))
                current =
                    (initial[ion] + reference[ion] - initial[electron] - reference[electron]) / 2
                @test maximum(x -> abs(Q(x) - current), st.J[component]) <= 2budget
            end
            @test all(s -> isequal(species[s].x[1], positions[s]), eachindex(species))
            @test all(ps -> all(iszero, ps.v[1]), species) && all(iszero, st.J[1])
            @test isequal(st.fields.B[1], bx)
            @test all(iszero, st.fields.E[1]) && all(a -> all(iszero, a), st.fields.B[2:3])
            @test all(iszero, st.charge_n) && all(iszero, st.charge_np1)
            @test r.gauss <= Q(256eps(T)) && r.continuity <= Q(256eps(T))
            @test r.energy_change == sum(r.defects)
            @test iszero(r.defects[2]) && iszero(r.curl_adjoint_work)
        end
    end
end
