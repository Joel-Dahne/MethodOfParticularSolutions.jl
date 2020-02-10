@testset "singular triangles" begin
    N = 16
    prec = 64

    RR = RealField(prec)

    old_prec = precision(BigFloat)
    setprecision(BigFloat, prec)

    results = RR.(["4.3 +/- 0.0813",
                   "5.16 +/- 7.74e-3",
                   "6.2 +/- 0.0849",
                   "6.78 +/- 8.61e-3"])

    for i in 1:4
        domain, u, interval = triangle(6 + i, RR)

        indices = filter(i -> numerator(domain.angles[i]) != 1, 1:3)

        λ, u = mps(domain, u, interval, N)

        @test overlaps(results[i], λ)
        @show λ
        @show Float64(radius(λ))
    end

    setprecision(BigFloat, old_prec)
end
