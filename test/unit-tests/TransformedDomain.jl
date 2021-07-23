@testset "TransformedDomain" begin
    has_rational_angles = MethodOfParticularSolutions.has_rational_angles
    angle = MethodOfParticularSolutions.angle
    angledivπ = MethodOfParticularSolutions.angledivπ
    anglesdivπ = MethodOfParticularSolutions.anglesdivπ
    vertices = MethodOfParticularSolutions.vertices
    orientation = MethodOfParticularSolutions.orientation
    center = MethodOfParticularSolutions.center

    parent = RealField(64)

    rotation = 1 // 3
    scaling = 2
    translation = [1, 2]
    triangle1 = Triangle(fmpq(1 // 3), fmpq(1 // 4), parent)
    domain1 = TransformedDomain(triangle1, rotation, scaling, translation)
    triangle2 = Triangle(parent(π) / 3, parent(π) / 4, parent)
    domain2 = TransformedDomain(triangle2, parent(π) * rotation, scaling, translation)

    for (triangle, domain) in [(triangle1, domain1), (triangle2, domain2)]
        @test has_rational_angles(domain) == ifelse(triangle isa Triangle{fmpq}, true, false)

        @test vertexindices(domain) == vertexindices(triangle)
        @test boundaries(domain) == boundaries(triangle)

        @test all(isequal(angle(domain, i), angles(domain)[i]) for i in boundaries(domain))
        @test all(
            isequal(angledivπ(domain, i), anglesdivπ(domain)[i]) for i in boundaries(domain)
        )
        @test isequal(angles(domain), angles(triangle))

        @test all(
            isequal(vertex(domain, i), vertices(domain)[i]) for i in boundaries(domain)
        )
        let (s, c) = sincospi(fmpq(rotation), parent)
            for i in boundaries(triangle)
                @test all(
                    overlaps.(
                        vertex(domain, i),
                        scaling * [c -s; s c] * vertex(triangle, i) + translation,
                    ),
                )
            end
        end

        if has_rational_angles(domain)
            @test isequal(orientation(domain, 1), orientation(triangle, 1) + rotation)
            @test isequal(orientation(domain, 2), orientation(triangle, 2) + rotation)
            @test isequal(orientation(domain, 3), orientation(triangle, 3) + rotation)
            @test isequal(
                orientation(domain, 1, reversed = true),
                orientation(triangle, 1, reversed = true) - rotation,
            )
            @test isequal(
                orientation(domain, 2, reversed = true),
                orientation(triangle, 2, reversed = true) - rotation,
            )
            @test isequal(
                orientation(domain, 3, reversed = true),
                orientation(triangle, 3, reversed = true) - rotation,
            )
        else
            @test overlaps(
                orientation(domain, 1),
                orientation(triangle, 1) + parent(π) * rotation,
            )
            @test overlaps(
                orientation(domain, 2),
                orientation(triangle, 2) + parent(π) * rotation,
            )
            @test overlaps(
                orientation(domain, 3),
                orientation(triangle, 3) + parent(π) * rotation,
            )
            @test overlaps(
                orientation(domain, 1, reversed = true),
                orientation(triangle, 1, reversed = true) - parent(π) * rotation,
            )
            @test overlaps(
                orientation(domain, 2, reversed = true),
                orientation(triangle, 2, reversed = true) - parent(π) * rotation,
            )
            @test overlaps(
                orientation(domain, 3, reversed = true),
                orientation(triangle, 3, reversed = true) - parent(π) * rotation,
            )
        end

        @test isequal(area(domain), scaling * area(triangle))

        let (s, c) = sincospi(fmpq(rotation), parent)
            @test all(
                overlaps.(
                    center(domain),
                    scaling * [c -s; s c] * center(triangle) + translation,
                ),
            )
        end
    end
end
