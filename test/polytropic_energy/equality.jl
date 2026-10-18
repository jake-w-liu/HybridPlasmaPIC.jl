module EqualityChecks
using Test, Random, HybridPlasmaPIC

const E = (
    prime = HybridPlasmaPIC._polytropic_sum_isprime,
    nextprime = HybridPlasmaPIC._polytropic_sum_nextprime,
    equal = HybridPlasmaPIC._polytropic_sum_equal,
)
const Q = Rational{BigInt}
println((VERSION,))
flush(stdout)
@testset "Deterministic modulus and exact multi-term controls" begin
    for n = 0:500
        expected = n >= 2 && all(n % d != 0 for d = 2:n-1)
        @test E.prime(BigInt(n)) == expected
    end
    modulus = BigInt(1)
    for _ = 1:50
        following = E.nextprime(modulus)
        @test following > modulus && E.prime(following)
        @test all(!E.prime(k) for k = modulus+1:following-1)
        modulus = following
    end
    rng = MersenneTwister(481253)
    for count = 2:7, p in (1, 2, 3, 7, 17), iteration = 1:3
        bases = Q[Q(rand(rng, 0:8), rand(rng, 1:9)) for _ = 1:count]
        target = sum(x -> x^p, bases)
        original = copy(bases)
        result = E.equal(bases, BigInt(p), target)
        @test result.equal
        @test !E.equal(bases, BigInt(p), target + Q(1, 97)).equal
        @test bases == original
        @test result.primes == 0 || result.certified_bits >= result.bound_bits
    end
    @test E.equal(Q[], BigInt(2), Q(0)).equal
    @test E.equal(Q[0, 0], BigInt(2), Q(0)).equal
    @test !E.equal(Q[0, 0], BigInt(2), Q(1)).equal
    @test !E.equal(Q[1, 2], BigInt(2), Q(-1)).equal
    @test_throws ArgumentError E.equal(Q[-1, 2], BigInt(2), Q(1))
    @test_throws ArgumentError E.equal(Q[1, 2], BigInt(0), Q(1))
end
flush(stdout)
@testset "Exact floating boundaries and huge exponent decisions" begin
    for T in (Float32, Float64)
        bits = precision(T)
        # Choose m=ceil(bits/2); multiplicity restores exactly half an ulp at one.
        m = cld(bits, 2)
        multiplicity = 1 << (2m - bits)
        bases = vcat(Q[1], fill(Q(1, BigInt(1) << m), multiplicity))
        boundary = (Q(one(T)) + Q(nextfloat(one(T)))) / 2
        result = E.equal(bases, BigInt(2), boundary)
        @test result.equal && result.primes > 0
        @test !E.equal(bases, BigInt(2), Q(one(T))).equal
        @test !E.equal(bases, BigInt(2), Q(nextfloat(one(T)))).equal
        println((T, result))
        flush(stdout)
    end
    huge = BigInt(10)^300
    for bases in (Q[1//2, 1//3], Q[1//7, 2//9, 3//10])
        @test all(x -> 0 < x < 1, bases) && sum(x -> x^2, bases) < 1 && huge > 2
        result = E.equal(bases, huge, Q(1))
        @test !result.equal
        @test result.max_modulus_bits < 16
        println((huge_exponent_bits = ndigits(huge; base = 2), result))
        flush(stdout)
    end
    @test E.equal(Q[1, 0, 1], huge, Q(2)).equal
    @test !E.equal(Q[2, 2], huge, Q(2)).equal
    bases = Q[3//5, 4//5]
    generator = (x for x in bases)
    @test E.equal(generator, BigInt(2), Q(1)).equal
end
flush(stdout)
@testset "Actual repeated-density resource measurements" begin
    for count in (32, 256, 4096)
        bases = [isodd(i) ? Q(3, 5) : Q(4, 5) for i = 1:count]
        target = Q(count, 2)
        E.equal(bases, BigInt(2), target)
        trial = @timed E.equal(bases, BigInt(2), target)
        @test trial.value.equal
        @test trial.value.certified_bits >= trial.value.bound_bits
        println((count, bytes = trial.bytes, seconds = trial.time, result = trial.value))
        flush(stdout)
    end
end
println("RATIONAL_POWER_SUM_EQUALITY_COMPLETE")
flush(stdout)

end
