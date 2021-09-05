@testset "Triangle" begin
    has_rational_angles = MethodOfParticularSolutions.has_rational_angles
    angle_raw = MethodOfParticularSolutions.angle_raw
    angle = MethodOfParticularSolutions.angle
    vertices = MethodOfParticularSolutions.vertices
    orientation_raw = MethodOfParticularSolutions.orientation_raw
    orientation = MethodOfParticularSolutions.orientation
    center = MethodOfParticularSolutions.center

    parent = RealField(64)

    domain1 = Triangle(fmpq(1 // 3), fmpq(1 // 4); parent)
    domain2 = Triangle(parent(π) / 3, parent(π) / 4; parent)
    domain3 = Triangle(1 // 3, 1 // 4)
    domain4 = Triangle(π / 3, π / 4)

    for domain in [domain1, domain2, domain3, domain4]
        @test typeof(typeof(domain)(domain)) == typeof(domain)

        @test has_rational_angles(domain) ==
              (domain isa Triangle{T,<:Union{Rational,fmpq}} where {T})

        @test vertexindices(domain) == 1:3
        @test boundaries(domain) == 1:3

        if has_rational_angles(domain)
            @test isequal(angle_raw(domain, 1), fmpq(1 // 3))
            @test isequal(angle_raw(domain, 2), fmpq(1 // 4))
            @test isequal(angle_raw(domain, 3), fmpq(5 // 12))
        else
            @test Float64(angle_raw(domain, 1)) ≈ π / 3
            @test Float64(angle_raw(domain, 2)) ≈ π / 4
            @test Float64(angle_raw(domain, 3)) ≈ 5π / 12
        end

        @test all(isequal(angle(domain, i), angles(domain)[i]) for i in boundaries(domain))
        @test all(isapprox.(Float64.(angles(domain)), [π / 3, π / 4, 5π / 12]))

        @test all(
            isequal(vertex(domain, i), vertices(domain)[i]) for i in boundaries(domain)
        )
        @test vertex(domain, 1) == [0, 0]
        @test vertex(domain, 2) == [1, 0]
        @test Float64.(vertex(domain, 3)) ≈
              sinpi(1 // 4) / sinpi(5 // 12) * [cospi(1 // 3), sinpi(1 // 3)]

        if has_rational_angles(domain)
            @test isequal(orientation_raw(domain, 1), fmpq(0))
            @test isequal(orientation_raw(domain, 2), fmpq(3 // 4))
            @test isequal(orientation_raw(domain, 3), fmpq(1 // 3 + 1))
            @test isequal(orientation_raw(domain, 1, reversed = true), fmpq(5 // 3))
            @test isequal(orientation_raw(domain, 2, reversed = true), fmpq(1))
            @test isequal(
                orientation_raw(domain, 3, reversed = true),
                fmpq(1 - 1 // 3 - 5 // 12),
            )
        else
            @test Float64(orientation_raw(domain, 1)) ≈ 0
            @test Float64(orientation_raw(domain, 2)) ≈ π * 3 // 4
            @test Float64(orientation_raw(domain, 3)) ≈ π * (1 // 3 + 1)
            @test Float64(orientation_raw(domain, 1, reversed = true)) ≈ π * 5 // 3
            @test Float64(orientation_raw(domain, 2, reversed = true)) ≈ π
            @test Float64(orientation_raw(domain, 3, reversed = true)) ≈
                  π * (1 - 1 // 3 - 5 // 12)
        end

        @test Float64(orientation(domain, 1)) ≈ 0
        @test Float64(orientation(domain, 2)) ≈ π * 3 // 4
        @test Float64(orientation(domain, 3)) ≈ π * (1 // 3 + 1)
        @test Float64(orientation(domain, 1, reversed = true)) ≈ π * 5 // 3
        @test Float64(orientation(domain, 2, reversed = true)) ≈ π
        @test Float64(orientation(domain, 3, reversed = true)) ≈ π * (1 - 1 // 3 - 5 // 12)

        @test Float64(area(domain)) ≈ sqrt(3) / 2(1 + sqrt(3))

        @test Float64.(center(domain)) ≈
              ([1, 0] + sinpi(1 // 4) / sinpi(5 // 12) * [cospi(1 // 3), sinpi(1 // 3)]) / 3

        # This only checks that the method runs and returns the
        # correct number of points, but not the actual values
        for i in vertexindices(domain)
            for distribution in [:linear, :chebyshev, :exponential, :root_exponential]
                points, boundary = boundary_points(domain, i, 10; distribution)
                @test length(points) == length(boundary) == 10
                @test boundary == fill(i, 10)
            end
        end

        points = interior_points(domain, 10)
        @test length(points) == 10
        @test all(∈(domain), points)
    end
end
