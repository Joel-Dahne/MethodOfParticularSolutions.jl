@testset "Polygon" begin
    has_rational_angles = MethodOfParticularSolutions.has_rational_angles
    angle_raw = MethodOfParticularSolutions.angle_raw
    angle = MethodOfParticularSolutions.angle
    vertices = MethodOfParticularSolutions.vertices
    orientation_raw = MethodOfParticularSolutions.orientation_raw
    orientation = MethodOfParticularSolutions.orientation
    center = MethodOfParticularSolutions.center

    parent = RealField(64)

    triangle1 = Triangle(fmpq(1 // 3), fmpq(1 // 4), parent)
    domain1 =
        Polygon([angle_raw(triangle1, i) for i = 1:3], collect(vertices(triangle1)), parent)
    triangle2 = Triangle(parent(π) / 3, parent(π) / 4, parent)
    domain2 = Polygon(collect(angles(triangle2)), collect(vertices(triangle2)), parent)

    for (triangle, domain) in [(triangle1, domain1), (triangle2, domain2)]
        @test has_rational_angles(domain) ==
              ifelse(triangle isa Triangle{fmpq}, true, false)

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
            # TODO: Not yet implemented for rational angles
            @test overlaps(orientation_raw(domain, 1), orientation_raw(triangle, 1))
            @test overlaps(orientation_raw(domain, 2), orientation_raw(triangle, 2))
            @test overlaps(orientation_raw(domain, 3), orientation_raw(triangle, 3))
            @test overlaps(
                orientation_raw(domain, 1, reversed = true),
                orientation_raw(triangle, 1, reversed = true),
            )
            @test overlaps(
                orientation_raw(domain, 2, reversed = true),
                orientation_raw(triangle, 2, reversed = true),
            )
            @test overlaps(
                orientation_raw(domain, 3, reversed = true),
                orientation_raw(triangle, 3, reversed = true),
            )
        end

        @test overlaps(orientation(domain, 1), orientation(triangle, 1))
        @test overlaps(orientation(domain, 2), orientation(triangle, 2))
        @test overlaps(orientation(domain, 3), orientation(triangle, 3))
        @test overlaps(
            orientation(domain, 1, reversed = true),
            orientation(triangle, 1, reversed = true),
        )
        @test overlaps(
            orientation(domain, 2, reversed = true),
            orientation(triangle, 2, reversed = true),
        )
        @test overlaps(
            orientation(domain, 3, reversed = true),
            orientation(triangle, 3, reversed = true),
        )

        @test overlaps(area(domain), area(triangle))

        @test all(overlaps.(center(domain), center(triangle)))
    end
end
