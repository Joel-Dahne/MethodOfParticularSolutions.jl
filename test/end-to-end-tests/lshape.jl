@testset "L-shape" begin
    Ns = 8:8:16
    optim_prec_final = 64

    parent = RealField(1024)
    result = parent("[9.6397238440219410527114592624 +/- 1.43e-26]")
    goalradius = 3e-2

    domain = LShape(RealField(128))
    u = LShapeEigenfunction(domain)
    interval = ball(parent(9.6), parent(0.5))

    println("Computing eigenvalue for the $domain")
    λs, _ = iteratemps(
        domain,
        u,
        interval,
        Ns,
        optim_prec_final = optim_prec_final,
        optim_prec_linear = true,
        show_trace = true,
    )
    rad = Float64(radius(λs[end]))

    @test overlaps(result, λs[end])
    @test rad < goalradius

    if rad < goalradius
        @printf "radius = %e < %e\n\n" rad goalradius
    else
        @printf "radius = %e >= %e\n\n" rad goalradius
    end
end
