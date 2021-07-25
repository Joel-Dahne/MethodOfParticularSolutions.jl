@testset "LinkedEigenfunction" begin
    parent = RealField(64)
    series = x -> arb_series(ArbPolyRing(parent, :x)(parent.(x)))
    triangle = Triangle(fmpq(1 // 3), fmpq(1 // 3), parent)

    u1 = LinkedEigenfunction(
        [StandaloneVertexEigenfunction{Float64,Rational{Int}}(triangle, i) for i = 1:3],
        excluded_boundaries = [BitSet([2, 3]), BitSet([1, 2]), BitSet([1, 3])],
    )
    u2 = LinkedEigenfunction(
        [StandaloneVertexEigenfunction{Float64,Float64}(triangle, i) for i = 1:3],
        excluded_boundaries = [BitSet([2, 3]), BitSet([1, 2]), BitSet([1, 3])],
    )
    u3 = LinkedEigenfunction(
        [StandaloneVertexEigenfunction{arb,fmpq}(triangle, i) for i = 1:3],
        excluded_boundaries = [BitSet([2, 3]), BitSet([1, 2]), BitSet([1, 3])],
    )
    u4 = LinkedEigenfunction(
        [StandaloneVertexEigenfunction{arb,arb}(triangle, i) for i = 1:3],
        excluded_boundaries = [BitSet([2, 3]), BitSet([1, 2]), BitSet([1, 3])],
    )

    for (u, S, T) in [
        (u1, Float64, Rational{Int}),
        (u2, Float64, Float64),
        (u3, arb, fmpq),
        (u4, arb, arb),
    ]
        @test u isa LinkedEigenfunction{S,T,StandaloneVertexEigenfunction{S,T}}
        @test length(u.us) == 3
        @test u.extra_coefficients == [1, 1, 1]
        @test u.excluded_boundaries == [BitSet([2, 3]), BitSet([1, 2]), BitSet([1, 3])]

        @test coefficients(set_eigenfunction!(u, [1, 2, 3])) == [1, 2, 3]
        @test all(i -> coefficients(u.us[i]) == [1, 2, 3], 1:3)

        if S == arb
            @test all(
                overlaps.(
                    u([1, 1], 1, 1:2)::Vector{S},
                    sum(v([1, 1], 1, 1:2) for v in u.us),
                ),
            )
        else
            @test u([1, 1], 1, 1:2)::Vector{S} == sum(v([1, 1], 1, 1:2) for v in u.us)
        end
    end

end
