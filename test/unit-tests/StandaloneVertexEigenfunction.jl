@testset "StandaloneVertexEigenfunction" begin
    has_rational_angles = MethodOfParticularSolutions.has_rational_angles
    angle_raw = MethodOfParticularSolutions.angle_raw
    orientation_raw = MethodOfParticularSolutions.orientation_raw
    nu = MethodOfParticularSolutions.nu
    coordinate_transformation = MethodOfParticularSolutions.coordinate_transformation

    parent = RealField(64)
    series = x -> arb_series(ArbPolyRing(parent, :x)(parent.(x)))

    u1 = StandaloneVertexEigenfunction([1.0, 2.0], 1 // 2, 2 // 3)
    u2 = StandaloneVertexEigenfunction{Float64,Rational{Int}}([1.0, 2.0], 1 // 2, 2 // 3)
    for u in [u1, u2]
        @test u isa StandaloneVertexEigenfunction{Float64,Rational{Int}}
        @test u.vertex == [1, 2]
        @test isequal(u.orientation, 1 // 2)
        @test isequal(u.θ, 2 // 3)
        @test u.stride == 1
        @test u.offset == 0
        @test u.reversed == false
        @test isempty(u.coefficients)
        @test isnothing(u.parent)

        @test coefficients(set_eigenfunction!(u, [1, 2, 3])) == [1, 2, 3]

        @test nu(u, 3)::Float64 == 9 // 2

        @test coordinate_transformation(u, [1, 3])::AbstractVector{Float64} == [1, 0]

        @test u([1, 1], 1, 1:2)::Vector{Float64} ≈ [-0.24029783912342725, 0]
    end

    u1 = StandaloneVertexEigenfunction([1.0, 2.0], π / 2, 2π / 3)
    u2 = StandaloneVertexEigenfunction{Float64,Float64}([1.0, 2.0], π / 2, 2π / 3)
    for u in [u1, u2]
        @test u isa StandaloneVertexEigenfunction{Float64,Float64}
        @test u.vertex == [1, 2]
        @test isequal(u.orientation, π / 2)
        @test isequal(u.θ, 2π / 3)
        @test u.stride == 1
        @test u.offset == 0
        @test u.reversed == false
        @test isempty(u.coefficients)
        @test isnothing(u.parent)

        @test coefficients(set_eigenfunction!(u, [1, 2, 3])) == [1, 2, 3]

        @test nu(u, 3)::Float64 ≈ 9 / 2

        @test coordinate_transformation(u, [1, 3])::AbstractVector{Float64} ≈ [1, 0]

        @test u([1, 1], 1, 1:2)::Vector{Float64} ≈ [-0.24029783912342725, 0]
    end

    u1 = StandaloneVertexEigenfunction(parent.([1.0, 2.0]), 1 // 2, 2 // 3)
    u2 = StandaloneVertexEigenfunction{arb,fmpq}(
        parent.([1.0, 2.0]),
        fmpq(1 // 2),
        fmpq(2 // 3),
    )
    for u in [u1, u2]
        @test u isa StandaloneVertexEigenfunction{arb,fmpq}
        @test u.vertex == [1, 2]
        @test isequal(u.orientation, fmpq(1 // 2))
        @test isequal(u.θ, fmpq(2 // 3))
        @test u.stride == 1
        @test u.offset == 0
        @test u.reversed == false
        @test isempty(u.coefficients)
        @test u.parent == parent

        @test coefficients(set_eigenfunction!(u, [1, 2, 3])) == [1, 2, 3]

        @test nu(u, 3)::arb == 9 // 2

        @test all(
            overlaps.(
                coordinate_transformation(u, [1, 3])::AbstractVector{arb},
                [parent(1), parent(0)],
            ),
        )

        @test Float64.(u([1, 1], 1, 1:2)::Vector{arb}) ≈ [-0.24029783912342725, 0]
        @test u([series([1]), series([1])], 1, 1:2) isa Vector{arb_series}
    end

    u1 = StandaloneVertexEigenfunction(parent.([1.0, 2.0]), parent(π) / 2, 2parent(π) / 3)
    u2 = StandaloneVertexEigenfunction{arb,arb}(
        parent.([1.0, 2.0]),
        parent(π) / 2,
        2parent(π) / 3,
    )
    for u in [u1, u2]
        @test u isa StandaloneVertexEigenfunction{arb,arb}
        @test u.vertex == [1, 2]
        @test isequal(u.orientation, parent(π) / 2)
        @test isequal(u.θ, 2parent(π) / 3)
        @test u.stride == 1
        @test u.offset == 0
        @test u.reversed == false
        @test isempty(u.coefficients)
        @test u.parent == parent

        @test coefficients(set_eigenfunction!(u, [1, 2, 3])) == [1, 2, 3]

        @test Float64(nu(u, 3)::arb) ≈ 9 / 2

        @test all(
            overlaps.(
                coordinate_transformation(u, [1, 3])::AbstractVector{arb},
                [parent(1), parent(0)],
            ),
        )

        @test Float64.(u([1, 1], 1, 1:2)::Vector{arb}) ≈ [-0.24029783912342725, 0]
        @test u([series([1]), series([1])], 1, 1:2) isa Vector{arb_series}
    end

    domain1 = Triangle(fmpq(1 // 3), fmpq(1 // 4), parent)
    domain2 = Triangle(parent(π) / 3, parent(π) / 4, parent)

    for domain in [domain1, domain2]
        u = StandaloneVertexEigenfunction(domain, 2)
        if has_rational_angles(domain)
            @test u isa StandaloneVertexEigenfunction{arb,fmpq}
        else
            @test u isa StandaloneVertexEigenfunction{arb,arb}
        end
        @test u.vertex == vertex(domain, 2)
        @test isequal(u.orientation, orientation_raw(domain, 2))
        @test isequal(u.θ, angle_raw(domain, 2))
        @test u.stride == 1
        @test u.offset == 0
        @test u.reversed == false
        @test isempty(u.coefficients)
        @test u.parent == domain.parent

        @test coefficients(set_eigenfunction!(u, [1, 2, 3])) == [1, 2, 3]

        if has_rational_angles(domain)
            @test nu(u, 3)::arb == 12
        else
            @test Float64(nu(u, 3)::arb) ≈ 12
        end

        @test Float64.(
            coordinate_transformation(u, vertex(domain, 3))::AbstractVector{arb},
        ) ≈ [cospi(1 / 6) / sinpi(5 / 12), 0]

        @test Float64.(u([0.1, 0.1], 1, 1:2)::Vector{arb}) ≈
              [0.0007197702171518307, -3.3141808172927214e-8]
        @test u([series([1]), series([1])], 1, 1:2) isa Vector{arb_series}

        if has_rational_angles(domain)
            u = StandaloneVertexEigenfunction{Float64,Rational{Int}}(domain, 2)
            @test u isa StandaloneVertexEigenfunction{Float64,Rational{Int}}
        else
            u = StandaloneVertexEigenfunction{Float64,Float64}(domain, 2)
            @test u isa StandaloneVertexEigenfunction{Float64,Float64}
        end
        @test u.vertex == Float64.(vertex(domain, 2))
        if has_rational_angles(domain)
            @test u.orientation == orientation_raw(domain, 2)
            @test u.θ == angle_raw(domain, 2)
        else
            @test u.orientation == Float64(orientation_raw(domain, 2))
            @test u.θ == Float64(angle_raw(domain, 2))
        end
        @test u.stride == 1
        @test u.offset == 0
        @test u.reversed == false
        @test isempty(u.coefficients)
        @test u.parent == domain.parent

        @test coefficients(set_eigenfunction!(u, [1, 2, 3])) == [1, 2, 3]

        if has_rational_angles(domain)
            @test nu(u, 3)::Float64 == 12
        else
            @test Float64(nu(u, 3)::Float64) ≈ 12
        end

        @test coordinate_transformation(
            u,
            Float64.(vertex(domain, 3)),
        )::AbstractVector{Float64} ≈ [cospi(1 / 6) / sinpi(5 / 12), 0]

        @test u([0.1, 0.1], 1, 1:2)::Vector{Float64} ≈
              [0.0007197702171518307, -3.3141808172927214e-8]
    end
end
