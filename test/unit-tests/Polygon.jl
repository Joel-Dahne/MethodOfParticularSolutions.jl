@testset "Polygon" begin
    has_rational_angles = MethodOfParticularSolutions.has_rational_angles
    angle_raw = MethodOfParticularSolutions.angle_raw
    angle = MethodOfParticularSolutions.angle
    vertices = MethodOfParticularSolutions.vertices
    orientation_raw = MethodOfParticularSolutions.orientation_raw
    orientation = MethodOfParticularSolutions.orientation
    center = MethodOfParticularSolutions.center

    parent = RealField(64)

    triangle1 = Triangle(fmpq(1 // 3), fmpq(1 // 4); parent)
    domain1 =
        Polygon([angle_raw(triangle1, i) for i = 1:3], collect(vertices(triangle1)); parent)
    triangle2 = Triangle(parent(π) / 3, parent(π) / 4; parent)
    domain2 = Polygon(collect(angles(triangle2)), collect(vertices(triangle2)); parent)
    triangle3 = Triangle(1 // 3, 1 // 4)
    domain3 = Polygon([angle_raw(triangle3, i) for i = 1:3], collect(vertices(triangle3)))
    triangle4 = Triangle(π / 3, π / 4)
    domain4 = Polygon(collect(angles(triangle4)), collect(vertices(triangle4)))

    for (triangle, domain, S) in [
        (triangle1, domain1, arb),
        (triangle2, domain2, arb),
        (triangle3, domain3, Float64),
        (triangle4, domain4, Float64),
    ]
        @test has_rational_angles(domain) ==
              (domain isa Polygon{S,<:Union{Rational,fmpq}} where {S})

        @test vertexindices(domain) == vertexindices(triangle)
        @test boundaries(domain) == boundaries(triangle)

        @test all(
            isequal(angle_raw(domain, i), angle_raw(triangle, i)) for
            i in vertexindices(domain)
        )

        @test all(isequal(angle(domain, i), angles(domain)[i]) for i in boundaries(domain))
        @test isequal(angles(domain), collect(angles(triangle)))

        @test all(
            isequal(vertex(domain, i), vertices(domain)[i]) for i in boundaries(domain)
        )
        @test isequal(vertices(domain), collect(vertices(triangle)))

        if !has_rational_angles(domain)
            # This method is only implemented for non-rational angles
            if S == arb
                for i in vertexindices(domain)
                    @test overlaps(orientation_raw(domain, i), orientation_raw(triangle, i))
                    @test overlaps(
                        orientation_raw(domain, i, reversed = true),
                        orientation_raw(triangle, i, reversed = true),
                    )
                end
            else
                for i in vertexindices(domain)
                    @test orientation_raw(domain, i) ≈ orientation_raw(triangle, i)
                    @test orientation_raw(domain, i, reversed = true) ≈
                          orientation_raw(triangle, i, reversed = true)
                end
            end
        end

        if S == arb
            for i in vertexindices(domain)
                @test overlaps(orientation(domain, i), orientation(triangle, i))
                @test overlaps(
                    orientation(domain, i, reversed = true),
                    orientation(triangle, i, reversed = true),
                )
            end

            @test overlaps(area(domain), area(triangle))

            @test all(overlaps.(center(domain), center(triangle)))
        else
            for i in vertexindices(domain)
                @test orientation(domain, i) ≈ orientation(triangle, i)
                @test orientation(domain, i, reversed = true) ≈
                      orientation(triangle, i, reversed = true)
            end

            @test area(domain) ≈ area(triangle)

            @test center(domain) ≈ center(triangle)
        end

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
