@testset "TransformedDomain" begin
    has_rational_angles = MethodOfParticularSolutions.has_rational_angles
    angle_raw = MethodOfParticularSolutions.angle_raw
    angle = MethodOfParticularSolutions.angle
    vertices = MethodOfParticularSolutions.vertices
    orientation_raw = MethodOfParticularSolutions.orientation_raw
    orientation = MethodOfParticularSolutions.orientation
    center = MethodOfParticularSolutions.center

    parent = RealField(64)

    rotation = 1 // 3
    scaling = 2
    translation = [1, 2]
    triangle1 = Triangle(fmpq(1 // 3), fmpq(1 // 4); parent)
    domain1 = TransformedDomain(triangle1, rotation, scaling, translation)
    triangle2 = Triangle(parent(π) / 3, parent(π) / 4; parent)
    domain2 = TransformedDomain(triangle2, parent(π) * rotation, scaling, translation)
    triangle3 = Triangle(1 // 3, 1 // 4)
    domain3 = TransformedDomain(triangle3, rotation, scaling, translation)
    triangle4 = Triangle(π / 3, π / 4)
    domain4 = TransformedDomain(triangle4, π * rotation, scaling, translation)

    for (triangle, domain) in [
        (triangle1, domain1),
        (triangle2, domain2),
        (triangle3, domain3),
        (triangle4, domain4),
    ]
        @test has_rational_angles(domain) ==
              (domain isa TransformedDomain{T,<:Union{Rational,fmpq}} where {T})

        @test vertexindices(domain) == vertexindices(triangle)
        @test boundaries(domain) == boundaries(triangle)

        @test all(
            isequal(angle_raw(domain, i), angle_raw(triangle, i)) for
            i in vertexindices(domain)
        )

        @test all(isequal(angle(domain, i), angles(domain)[i]) for i in boundaries(domain))
        @test isequal(angles(domain), angles(triangle))

        @test all(
            isequal(vertex(domain, i), vertices(domain)[i]) for i in boundaries(domain)
        )
        if typeof(domain) <: TransformedDomain{arb}
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
        else
            let (s, c) = sincospi(rotation)
                for i in boundaries(triangle)
                    @test vertex(domain, i) ≈
                          scaling * [c -s; s c] * vertex(triangle, i) + translation
                end
            end
        end

        if has_rational_angles(domain)
            @test orientation_raw(domain, 1) == orientation_raw(triangle, 1) + rotation
            @test orientation_raw(domain, 2) == orientation_raw(triangle, 2) + rotation
            @test orientation_raw(domain, 3) == orientation_raw(triangle, 3) + rotation

            @test orientation_raw(domain, 1, reversed = true) ==
                  orientation_raw(triangle, 1, reversed = true) - rotation

            @test orientation_raw(domain, 2, reversed = true) ==
                  orientation_raw(triangle, 2, reversed = true) - rotation

            @test orientation_raw(domain, 3, reversed = true) ==
                  orientation_raw(triangle, 3, reversed = true) - rotation
        elseif typeof(domain) <: TransformedDomain{arb}
            @test overlaps(
                orientation_raw(domain, 1),
                orientation_raw(triangle, 1) + parent(π) * rotation,
            )
            @test overlaps(
                orientation_raw(domain, 2),
                orientation_raw(triangle, 2) + parent(π) * rotation,
            )
            @test overlaps(
                orientation_raw(domain, 3),
                orientation_raw(triangle, 3) + parent(π) * rotation,
            )
            @test overlaps(
                orientation_raw(domain, 1, reversed = true),
                orientation_raw(triangle, 1, reversed = true) - parent(π) * rotation,
            )
            @test overlaps(
                orientation_raw(domain, 2, reversed = true),
                orientation_raw(triangle, 2, reversed = true) - parent(π) * rotation,
            )
            @test overlaps(
                orientation_raw(domain, 3, reversed = true),
                orientation_raw(triangle, 3, reversed = true) - parent(π) * rotation,
            )
        else
            @test orientation_raw(domain, 1) == orientation_raw(triangle, 1) + π * rotation
            @test orientation_raw(domain, 2) == orientation_raw(triangle, 2) + π * rotation
            @test orientation_raw(domain, 3) == orientation_raw(triangle, 3) + π * rotation

            @test orientation_raw(domain, 1, reversed = true) ==
                  orientation_raw(triangle, 1, reversed = true) - π * rotation

            @test orientation_raw(domain, 2, reversed = true) ==
                  orientation_raw(triangle, 2, reversed = true) - π * rotation

            @test orientation_raw(domain, 3, reversed = true) ==
                  orientation_raw(triangle, 3, reversed = true) - π * rotation
        end

        @test isequal(area(domain), scaling * area(triangle))

        if typeof(domain) <: TransformedDomain{arb}
            let (s, c) = sincospi(fmpq(rotation), parent)
                @test all(
                    overlaps.(
                        center(domain),
                        scaling * [c -s; s c] * center(triangle) + translation,
                    ),
                )
            end
        else
            let (s, c) = sincospi(rotation)
                @test center(domain) ≈
                      scaling * [c -s; s c] * center(triangle) + translation
            end
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
