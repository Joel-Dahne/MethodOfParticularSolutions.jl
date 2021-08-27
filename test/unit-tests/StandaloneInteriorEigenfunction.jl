@testset "StandaloneInteriorEigenfunction" begin
    has_rational_angles = MethodOfParticularSolutions.has_rational_angles
    angle_raw = MethodOfParticularSolutions.angle_raw
    orientation_raw = MethodOfParticularSolutions.orientation_raw
    coordinate_transformation = MethodOfParticularSolutions.coordinate_transformation

    parent = RealField(64)
    series = x -> arb_series(ArbPolyRing(parent, :x)(parent.(x)))

    u1 = StandaloneInteriorEigenfunction([1.0, 2.0], 1 // 2)
    u2 = StandaloneInteriorEigenfunction{Float64,Rational{Int}}([1.0, 2.0], 1 // 2)
    for u in [u1, u2]
        @test u isa StandaloneInteriorEigenfunction{Float64,Rational{Int}}
        @test u.vertex == [1, 2]
        @test isequal(u.orientation, 1 // 2)
        @test u.stride == 1
        @test u.offset == 0
        @test !u.even && !u.odd
        @test isempty(u.coefficients)
        @test isnothing(u.parent)

        @test coefficients(set_eigenfunction!(u, [1, 2, 3])) == [1, 2, 3]

        @test coordinate_transformation(u, [1, 3])::AbstractVector{Float64} == [1, 0]

        @test u([1, 1], 1, 1:2)::Vector{Float64} ≈ [0.7651976865579666, 0]
    end

    u1 = StandaloneInteriorEigenfunction([1.0, 2.0], π / 2)
    u2 = StandaloneInteriorEigenfunction{Float64,Float64}([1.0, 2.0], π / 2)
    for u in [u1, u2]
        @test u isa StandaloneInteriorEigenfunction{Float64,Float64}
        @test u.vertex == [1, 2]
        @test isequal(u.orientation, π / 2)
        @test u.stride == 1
        @test u.offset == 0
        @test !u.even && !u.odd
        @test isempty(u.coefficients)
        @test isnothing(u.parent)

        @test coefficients(set_eigenfunction!(u, [1, 2, 3])) == [1, 2, 3]

        @test coordinate_transformation(u, [1, 3])::AbstractVector{Float64} ≈ [1, 0]

        @test u([1, 1], 1, 1:2)::Vector{Float64} ≈ [0.7651976865579666, 0]
    end

    u1 = StandaloneInteriorEigenfunction(parent.([1.0, 2.0]), 1 // 2)
    u2 = StandaloneInteriorEigenfunction{arb,fmpq}(parent.([1.0, 2.0]), fmpq(1 // 2))
    for u in [u1, u2]
        @test u isa StandaloneInteriorEigenfunction{arb,fmpq}
        @test u.vertex == [1, 2]
        @test isequal(u.orientation, fmpq(1 // 2))
        @test u.stride == 1
        @test u.offset == 0
        @test !u.even && !u.odd
        @test isempty(u.coefficients)
        @test u.parent == parent

        @test coefficients(set_eigenfunction!(u, [1, 2, 3])) == [1, 2, 3]

        @test all(
            overlaps.(
                coordinate_transformation(u, [1, 3])::AbstractVector{arb},
                [parent(1), parent(0)],
            ),
        )

        @test Float64.(u([1, 1], 1, 1:2)::Vector{arb}) ≈ [0.7651976865579666, 0]
        @test u([series([1]), series([1])], 1, 1:2) isa Vector{arb_series}
    end

    u1 = StandaloneInteriorEigenfunction(parent.([1.0, 2.0]), parent(π) / 2)
    u2 = StandaloneInteriorEigenfunction{arb,arb}(parent.([1.0, 2.0]), parent(π) / 2)
    for u in [u1, u2]
        @test u isa StandaloneInteriorEigenfunction{arb,arb}
        @test u.vertex == [1, 2]
        @test isequal(u.orientation, parent(π) / 2)
        @test u.stride == 1
        @test u.offset == 0
        @test !u.even && !u.odd
        @test isempty(u.coefficients)
        @test u.parent == parent

        @test coefficients(set_eigenfunction!(u, [1, 2, 3])) == [1, 2, 3]

        @test all(
            overlaps.(
                coordinate_transformation(u, [1, 3])::AbstractVector{arb},
                [parent(1), parent(0)],
            ),
        )

        @test Float64.(u([1, 1], 1, 1:2)::Vector{arb}) ≈ [0.7651976865579666, 0]
        @test u([series([1]), series([1])], 1, 1:2) isa Vector{arb_series}
    end

    domain1 = Triangle(fmpq(1 // 3), fmpq(1 // 4); parent)
    domain2 = Triangle(parent(π) / 3, parent(π) / 4; parent)

    for domain in [domain1, domain2]
        u = StandaloneInteriorEigenfunction(domain)
        if has_rational_angles(domain)
            @test u isa StandaloneInteriorEigenfunction{arb,fmpq}
        else
            @test u isa StandaloneInteriorEigenfunction{arb,arb}
        end
        @test Float64.(u.vertex) ≈ Float64.(center(domain))
        @test iszero(u.orientation)
        @test u.stride == 1
        @test u.offset == 0
        @test !u.even && !u.odd
        @test isempty(u.coefficients)
        @test u.parent == domain.parent

        @test coefficients(set_eigenfunction!(u, [1, 2, 3])) == [1, 2, 3]

        @test Float64.(coordinate_transformation(u, center(domain))::AbstractVector{arb},) ≈
              [0, 0]

        @test Float64.(u([0.1, 0.1], 1, 1:2)::Vector{arb}) ≈
              [0.9656340100141073, -0.0547032144426457]
        @test u([series([1]), series([1])], 1, 1:2) isa Vector{arb_series}

        if has_rational_angles(domain)
            u = StandaloneInteriorEigenfunction{Float64,Rational{Int}}(domain)
            @test u isa StandaloneInteriorEigenfunction{Float64,Rational{Int}}
        else
            u = StandaloneInteriorEigenfunction{Float64,Float64}(domain)
            @test u isa StandaloneInteriorEigenfunction{Float64,Float64}
        end
        @test u.vertex ≈ Float64.(center(domain))
        @test iszero(u.orientation)
        @test u.stride == 1
        @test u.offset == 0
        @test !u.even && !u.odd
        @test isempty(u.coefficients)
        @test u.parent == domain.parent

        @test coefficients(set_eigenfunction!(u, [1, 2, 3])) == [1, 2, 3]

        @test coordinate_transformation(
            u,
            Float64.(center(domain)),
        )::AbstractVector{Float64} ≈ [0, 0]

        @test u([0.1, 0.1], 1, 1:2)::Vector{Float64} ≈
              [0.9656340100141073, -0.0547032144426457]
    end
end
