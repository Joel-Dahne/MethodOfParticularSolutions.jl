@testset "StandaloneLightningEigenfunction" begin
    has_rational_angles = MethodOfParticularSolutions.has_rational_angles
    angle_raw = MethodOfParticularSolutions.angle_raw
    orientation_raw = MethodOfParticularSolutions.orientation_raw

    parent = RealField(64)

    u1 = StandaloneLightningEigenfunction([1.0, 2.0], 1 // 2, 2 // 3)
    u2 = StandaloneLightningEigenfunction{Float64,Rational{Int}}([1.0, 2.0], 1 // 2, 2 // 3)
    for u in [u1, u2]
        @test u isa StandaloneLightningEigenfunction{Float64,Rational{Int}}
        @test u.vertex == [1, 2]
        @test isequal(u.orientation, 1 // 2)
        @test isequal(u.θ, 2 // 3)
        @test u.l == 1
        @test u.σ == 4
        @test u.even == false
        @test u.odd == false
        @test u.reversed == false
        @test isempty(u.coefficients)
        @test isnothing(u.parent)
    end

    u1 = StandaloneLightningEigenfunction([1.0, 2.0], π / 2, 2π / 3)
    u2 = StandaloneLightningEigenfunction{Float64,Float64}([1.0, 2.0], π / 2, 2π / 3)
    for u in [u1, u2]
        @test u isa StandaloneLightningEigenfunction{Float64,Float64}
        @test u.vertex == [1, 2]
        @test isequal(u.orientation, π / 2)
        @test isequal(u.θ, 2π / 3)
        @test u.l == 1
        @test u.σ == 4
        @test u.even == false
        @test u.odd == false
        @test u.reversed == false
        @test isempty(u.coefficients)
        @test isnothing(u.parent)
    end

    u1 = StandaloneLightningEigenfunction(parent.([1.0, 2.0]), 1 // 2, 2 // 3)
    u2 = StandaloneLightningEigenfunction{arb,fmpq}(
        parent.([1.0, 2.0]),
        fmpq(1 // 2),
        fmpq(2 // 3),
    )
    for u in [u1, u2]
        @test u isa StandaloneLightningEigenfunction{arb,fmpq}
        @test u.vertex == [1, 2]
        @test isequal(u.orientation, fmpq(1 // 2))
        @test isequal(u.θ, fmpq(2 // 3))
        @test u.l == 1
        @test u.σ == 4
        @test u.even == false
        @test u.odd == false
        @test u.reversed == false
        @test isempty(u.coefficients)
        @test u.parent == parent
    end

    u1 = StandaloneLightningEigenfunction(
        parent.([1.0, 2.0]),
        parent(1 // 2),
        parent(2 // 3),
    )
    u2 = StandaloneLightningEigenfunction{arb,arb}(
        parent.([1.0, 2.0]),
        parent(1 // 2),
        parent(2 // 3),
    )
    for u in [u1, u2]
        @test u isa StandaloneLightningEigenfunction{arb,arb}
        @test u.vertex == [1, 2]
        @test isequal(u.orientation, parent(1 // 2))
        @test isequal(u.θ, parent(2 // 3))
        @test u.l == 1
        @test u.σ == 4
        @test u.even == false
        @test u.odd == false
        @test u.reversed == false
        @test isempty(u.coefficients)
        @test u.parent == parent
    end

    domain1 = Triangle(fmpq(1 // 3), fmpq(1 // 4), parent)
    domain2 = Triangle(parent(π) / 3, parent(π) / 4, parent)

    for domain in [domain1, domain2]
        u = StandaloneLightningEigenfunction(domain, 2)
        if has_rational_angles(domain)
            @test u isa StandaloneLightningEigenfunction{arb,fmpq}
        else
            @test u isa StandaloneLightningEigenfunction{arb,arb}
        end
        @test u.vertex == vertex(domain, 2)
        @test isequal(u.orientation, orientation_raw(domain, 2))
        @test isequal(u.θ, angle_raw(domain, 2))
        @test u.l == 1
        @test u.σ == 4
        @test u.even == false
        @test u.odd == false
        @test u.reversed == false
        @test isempty(u.coefficients)
        @test u.parent == domain.parent

        if has_rational_angles(domain)
            u = StandaloneLightningEigenfunction{Float64,Rational{Int}}(domain, 2)
            @test u isa StandaloneLightningEigenfunction{Float64,Rational{Int}}
        else
            u = StandaloneLightningEigenfunction{Float64,Float64}(domain, 2)
            @test u isa StandaloneLightningEigenfunction{Float64,Float64}
        end
        @test u.vertex == Float64.(vertex(domain, 2))
        if has_rational_angles(domain)
            @test u.orientation == orientation_raw(domain, 2)
            @test u.θ == angle_raw(domain, 2)
        else
            @test u.orientation == Float64(orientation_raw(domain, 2))
            @test u.θ == Float64(angle_raw(domain, 2))
        end
        @test u.l == 1
        @test u.σ == 4
        @test u.even == false
        @test u.odd == false
        @test u.reversed == false
        @test isempty(u.coefficients)
        @test u.parent == domain.parent
    end
end
