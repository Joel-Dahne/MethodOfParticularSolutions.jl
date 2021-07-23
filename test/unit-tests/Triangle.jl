@testset "Triangle" begin
    has_rational_angles = MethodOfParticularSolutions.has_rational_angles
    angle = MethodOfParticularSolutions.angle
    angledivπ = MethodOfParticularSolutions.angledivπ
    anglesdivπ = MethodOfParticularSolutions.anglesdivπ
    vertices = MethodOfParticularSolutions.vertices
    orientation = MethodOfParticularSolutions.orientation
    center = MethodOfParticularSolutions.center

    parent = RealField(64)

    domain1 = Triangle(fmpq(1 // 3), fmpq(1 // 4), parent)
    domain2 = Triangle(parent(π) / 3, parent(π) / 4, parent)

    for domain in [domain1, domain2]
        @test has_rational_angles(domain) == ifelse(domain isa Triangle{fmpq}, true, false)

        @test vertexindices(domain) == 1:3
        @test boundaries(domain) == 1:3

        @test all(isequal(angle(domain, i), angles(domain)[i]) for i in boundaries(domain))
        @test all(
            isequal(angledivπ(domain, i), anglesdivπ(domain)[i]) for i in boundaries(domain)
        )
        @test all(isapprox.(Float64.(angles(domain)), [π / 3, π / 4, 5π / 12]))

        if has_rational_angles(domain)
            @test all(isequal.(anglesdivπ(domain), fmpq.([1 // 3, 1 // 4, 5 // 12])))
        else
            @test all(isapprox.(Float64.(anglesdivπ(domain)), [1 // 3, 1 // 4, 5 // 12]))
        end

        @test all(
            isequal(vertex(domain, i), vertices(domain)[i]) for i in boundaries(domain)
        )
        @test vertex(domain, 1) == [0, 0]
        @test vertex(domain, 2) == [1, 0]
        @test Float64.(vertex(domain, 3)) ≈
              sinpi(1 // 4) / sinpi(5 // 12) * [cospi(1 // 3), sinpi(1 // 3)]

        if has_rational_angles(domain)
            @test isequal(orientation(domain, 1), fmpq(0))
            @test isequal(orientation(domain, 2), fmpq(3 // 4))
            @test isequal(orientation(domain, 3), fmpq(1 // 3 + 1))
            @test isequal(orientation(domain, 1, reversed = true), fmpq(5 // 3))
            @test isequal(orientation(domain, 2, reversed = true), fmpq(1))
            @test isequal(
                orientation(domain, 3, reversed = true),
                fmpq(1 - 1 // 3 - 5 // 12),
            )
        else
            @test Float64(orientation(domain, 1)) ≈ 0
            @test Float64(orientation(domain, 2)) ≈ π * 3 // 4
            @test Float64(orientation(domain, 3)) ≈ π * (1 // 3 + 1)
            @test Float64(orientation(domain, 1, reversed = true)) ≈ π * 5 // 3
            @test Float64(orientation(domain, 2, reversed = true)) ≈ π
            @test Float64(orientation(domain, 3, reversed = true)) ≈
                  π * (1 - 1 // 3 - 5 // 12)
        end

        @test Float64(area(domain)) ≈ sqrt(3) / 2(1 + sqrt(3))

        @test Float64.(center(domain)) ≈
              ([1, 0] + sinpi(1 // 4) / sinpi(5 // 12) * [cospi(1 // 3), sinpi(1 // 3)]) / 3
    end
end
