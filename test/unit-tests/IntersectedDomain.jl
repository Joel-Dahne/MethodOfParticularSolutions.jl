@testset "IntersectedDomain" begin
    has_rational_angles = MethodOfParticularSolutions.has_rational_angles
    angle_raw = MethodOfParticularSolutions.angle_raw
    angle = MethodOfParticularSolutions.angle
    vertices = MethodOfParticularSolutions.vertices
    orientation_raw = MethodOfParticularSolutions.orientation_raw
    orientation = MethodOfParticularSolutions.orientation
    center = MethodOfParticularSolutions.center

    parent = RealField(64)

    exterior1 = Triangle(fmpq(1 // 3), fmpq(1 // 4); parent)
    interior1 = Polygon(
        [1 // 4, 1 // 4, 1 // 4, 1 // 4],
        [[0.25, 0.1], [0.5, 0.1], [0.5, 0.2], [0.25, 0.2]];
        parent,
    )
    domain1 = IntersectedDomain(exterior1, [interior1])
    exterior2 = Triangle(parent(π) / 3, parent(π) / 4; parent)
    interior2 = Polygon(
        [parent(π) / 4, parent(π) / 4, parent(π) / 4, parent(π) / 4],
        [[0.25, 0.1], [0.5, 0.1], [0.5, 0.2], [0.25, 0.2]];
        parent,
    )
    domain2 = IntersectedDomain(exterior2, [interior2])

    for (exterior, interior, domain) in
        [(exterior1, interior1, domain1), (exterior2, interior2, domain2)]
        @test has_rational_angles(domain) ==
              ifelse(exterior isa Triangle{fmpq}, true, false)

        @test vertexindices(domain) == 1:7
        @test boundaries(domain) == 1:7

        @test all(isequal(angle_raw(domain, i), angle_raw(exterior, i)) for i = 1:3)
        @test all(isequal(angle_raw(domain, i + 3), angle_raw(interior, i)) for i = 1:4)

        @test all(isequal(angle(domain, i), angles(domain)[i]) for i in boundaries(domain))
        @test all(isequal(angle(domain, i), angle(exterior, i)) for i = 1:3)
        @test all(isequal(angle(domain, i + 3), angle(interior, i)) for i = 1:4)

        @test all(
            isequal(vertex(domain, i), vertices(domain)[i]) for i in boundaries(domain)
        )
        @test all(isequal(vertex(domain, i), vertex(exterior, i)) for i = 1:3)
        @test all(isequal(vertex(domain, i + 3), vertex(interior, i)) for i = 1:4)

        @test all(
            isequal(orientation_raw(domain, i), orientation_raw(exterior, i)) for i = 1:3
        )
        # TODO: Not implemented for rational polygons
        #@test all(isequal(orientation_raw(domain, i + 3), orientation_raw(interior, i)) for i = 1:4)
        @test all(
            isequal(
                orientation_raw(domain, i, reversed = true),
                orientation_raw(exterior, i, reversed = true),
            ) for i = 1:3
        )
        # TODO: Not implemented for rational polygons
        #@test all(
        #    isequal(
        #        orientation_raw(domain, i + 3, reversed = true),
        #        orientation_raw(interior, i, reversed = true),
        #    ) for i = 1:4
        #)

        @test all(isequal(orientation(domain, i), orientation(exterior, i)) for i = 1:3)
        # TODO: Not implemented for rational polygons
        #@test all(isequal(orientation(domain, i + 3), orientation(interior, i)) for i = 1:4)
        @test all(
            isequal(
                orientation(domain, i, reversed = true),
                orientation(exterior, i, reversed = true),
            ) for i = 1:3
        )
        # TODO: Not implemented for rational polygons
        #@test all(
        #    isequal(
        #        orientation(domain, i + 3, reversed = true),
        #        orientation(interior, i, reversed = true),
        #    ) for i = 1:4
        #)

        @test isequal(area(domain), area(exterior) - area(interior))

        @test isequal(center(domain), center(exterior))
    end
end
