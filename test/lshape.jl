@testset "L-shape" begin
    N = 16
    RR = RealField(64)

    result = RR("9.6397238440219410527114592624 +/- 1.43e-26")

    domain = LShape(RR)
    u = LShapeEigenfunction(domain)
    interval = ball(RR(9.6), RR(0.5))

    λ, u = setprecision(BigFloat, prec(RR)) do
        mps(domain, u, interval, N)
    end

    @test overlaps(result, λ)
    @show λ
    @show Float64(radius(λ))
end
