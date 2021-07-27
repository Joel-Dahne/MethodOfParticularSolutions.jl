@testset "CombinedEigenfunction" begin
    basis_function = MethodOfParticularSolutions.basis_function

    parent = RealField(64)
    series = x -> arb_series(ArbPolyRing(parent, :x)(parent.(x)))
    triangle = Triangle(fmpq(1 // 5), fmpq(1 // 5), parent)

    u1 = CombinedEigenfunction{Float64,Rational{Int}}(
        triangle,
        [StandaloneVertexEigenfunction{Float64,Rational{Int}}(triangle, i) for i = 1:3],
        [1, 1, 2],
        BitSet([3]),
        Dict(i => BitSet([i]) for i = 1:3),
    )
    u2 = CombinedEigenfunction{Float64,Float64}(
        triangle,
        [StandaloneVertexEigenfunction{Float64,Float64}(triangle, i) for i = 1:3],
        [1, 1, 2],
        BitSet([3]),
        Dict(i => BitSet([i]) for i = 1:3),
    )
    u3 = CombinedEigenfunction{arb,fmpq}(
        triangle,
        [StandaloneVertexEigenfunction{arb,fmpq}(triangle, i) for i = 1:3],
        [1, 1, 2],
        BitSet([3]),
        Dict(i => BitSet([i]) for i = 1:3),
    )
    u4 = CombinedEigenfunction{arb,arb}(
        triangle,
        [StandaloneVertexEigenfunction{arb,arb}(triangle, i) for i = 1:3],
        [1, 1, 2],
        BitSet([3]),
        Dict(i => BitSet([i]) for i = 1:3),
    )

    for (u, S, T) in [
        (u1, Float64, Rational{Int}),
        (u2, Float64, Float64),
        (u3, arb, fmpq),
        (u4, arb, arb),
    ]
        @test u isa CombinedEigenfunction{S,T}
        @test length(u.us) == 3
        @test u.orders == [1, 1, 2]
        @test u.even_boundaries == BitSet([3])
        @test u.boundary_to_us == Dict(i => BitSet([i]) for i = 1:3)

        @test coefficients(set_eigenfunction!(u, [1, 2, 3, 4, 5, 6, 7])) ==
              [1, 2, 3, 4, 5, 6, 7]
        @test coefficients(u.us[1]) == [1, 5]
        @test coefficients(u.us[2]) == [2, 6]
        @test coefficients(u.us[3]) == [3, 4, 7]

        @test basis_function(u, 1) == (1, 1)
        @test basis_function(u, 2) == (2, 1)
        @test basis_function(u, 3) == (3, 1)
        @test basis_function(u, 4) == (3, 2)
        @test basis_function(u, 5) == (1, 2)
        @test basis_function(u, 1:5) == [1:2, 1:1, 1:2]
        @test basis_function(u, 1:10) == [1:3, 1:3, 1:4]
        @test basis_function(u, 5:10) == [2:3, 2:3, 3:4]
        @test basis_function(u, 10:10) == [3:2, 3:3, 4:3]

        if S == arb
            @test overlaps(u([1, 1], 1)::S, sum(v([1, 1], 1) for v in u.us))
            @test isequal(u([1, 1], 1, 2)::S, u.us[2]([1, 1], 1, 1))
            @test all(
                isequal.(
                    u([1, 1], 1, 1:5)::Vector{S},
                    [
                        u.us[1]([1, 1], 1, 1),
                        u.us[2]([1, 1], 1, 1),
                        u.us[3]([1, 1], 1, 1),
                        u.us[3]([1, 1], 1, 2),
                        u.us[1]([1, 1], 1, 2),
                    ],
                ),
            )

            @test overlaps(
                (u(series.([1, 1]), 1)::arb_series)[0],
                sum(v(series.([1, 1]), 1) for v in u.us)[0],
            )
            @test isequal(
                (u(series.([1, 1]), 1, 2)::arb_series)[0],
                u.us[2](series.([1, 1]), 1, 1)[0],
            )
            @test all(
                isequal.(
                    map(
                        x -> getindex(x, 0),
                        u(series.([1, 1]), 1, 1:5)::Vector{arb_series},
                    ),
                    [
                        u.us[1](series.([1, 1]), 1, 1)[0],
                        u.us[2](series.([1, 1]), 1, 1)[0],
                        u.us[3](series.([1, 1]), 1, 1)[0],
                        u.us[3](series.([1, 1]), 1, 2)[0],
                        u.us[1](series.([1, 1]), 1, 2)[0],
                    ],
                ),
            )
        else
            @test u([1, 1], 1)::S == sum(v([1, 1], 1) for v in u.us)
            @test u([1, 1], 1, 2)::S == u.us[2]([1, 1], 1, 1)
            @test u([1, 1], 1, 1:5)::Vector{S} == [
                u.us[1]([1, 1], 1, 1),
                u.us[2]([1, 1], 1, 1),
                u.us[3]([1, 1], 1, 1),
                u.us[3]([1, 1], 1, 2),
                u.us[1]([1, 1], 1, 2),
            ]
        end
    end

end
