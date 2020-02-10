begin
    N = 16
    prec = 64

    RR = RealField(64)

    old_prec = precision(BigFloat)
    setprecision(BigFloat, prec)

    result = RR("9.6397238440219410527114592624 +/- 1.43e-26")

    domain = LShape(RR)
    u = LShapeEigenfunction(domain)
    interval = ball(RR(9.6), RR(0.5))

    λ, u = mps(domain, u, interval, N)

    @test overlaps(result, λ)
    @show λ
    @show Float64(radius(λ))

    setprecision(BigFloat, old_prec)
end
