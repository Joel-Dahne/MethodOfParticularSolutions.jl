@testset "Polygon" begin
    angle = MethodOfParticularSolutions.angle
    angledivπ = MethodOfParticularSolutions.angledivπ
    anglesdivπ = MethodOfParticularSolutions.anglesdivπ
    vertices = MethodOfParticularSolutions.vertices
    orientation = MethodOfParticularSolutions.orientation
    center = MethodOfParticularSolutions.center

    parent = RealField(64)

    triangle1 = Triangle(fmpq(1 // 3), fmpq(1 // 4), parent)
    domain1 = Polygon(collect(anglesdivπ(triangle1)), collect(vertices(triangle1)), parent)
    triangle2 = Triangle(parent(π) / 3, parent(π) / 4, parent)
    domain2 = Polygon(collect(angles(triangle2)), collect(vertices(triangle2)), parent)

    for (triangle, domain) in [(triangle1, domain1), (triangle2, domain2)]
        @test boundaries(domain) == boundaries(triangle)

        @test all(isequal(angle(domain, i), angles(domain)[i]) for i in boundaries(domain))
        @test all(
            isequal(angledivπ(domain, i), anglesdivπ(domain)[i]) for i in boundaries(domain)
        )
        @test isequal(angles(domain), collect(angles(triangle)))

        @test all(
            isequal(vertex(domain, i), vertices(domain)[i]) for i in boundaries(domain)
        )
        @test isequal(vertices(domain), collect(vertices(triangle)))

        if triangle isa Triangle{fmpq}
            @test overlaps(
                orientation(domain, 1),
                parent(π) * parent(orientation(triangle, 1)),
            )
            @test overlaps(
                orientation(domain, 2),
                parent(π) * parent(orientation(triangle, 2)),
            )
            @test overlaps(
                orientation(domain, 3),
                parent(π) * parent(orientation(triangle, 3)),
            )
            @test overlaps(
                orientation(domain, 1, reversed = true),
                parent(π) * parent(orientation(triangle, 1, reversed = true)),
            )
            @test overlaps(
                orientation(domain, 2, reversed = true),
                parent(π) * parent(orientation(triangle, 2, reversed = true)),
            )
            @test overlaps(
                orientation(domain, 3, reversed = true),
                parent(π) * parent(orientation(triangle, 3, reversed = true)),
            )
        else
            @test overlaps(orientation(domain, 1), parent(orientation(triangle, 1)))
            @test overlaps(orientation(domain, 2), parent(orientation(triangle, 2)))
            @test overlaps(orientation(domain, 3), parent(orientation(triangle, 3)))
            @test overlaps(
                orientation(domain, 1, reversed = true),
                parent(orientation(triangle, 1, reversed = true)),
            )
            @test overlaps(
                orientation(domain, 2, reversed = true),
                parent(orientation(triangle, 2, reversed = true)),
            )
            @test overlaps(
                orientation(domain, 3, reversed = true),
                parent(orientation(triangle, 3, reversed = true)),
            )
        end

        @test overlaps(area(domain), area(triangle))

        @test all(overlaps.(center(domain), center(triangle)))
    end
end
