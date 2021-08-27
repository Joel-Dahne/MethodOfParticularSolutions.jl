@testset "triangle" begin
    println(
        "Computing eigenvalues for the equilateral triangle with different types of eigenfunctions",
    )

    # Take an equilateral triangle
    triangle = Triangle(fmpq(1 // 3), fmpq(1 // 3), parent = RealField(64))

    # Take two different kind of eigenfunctions
    u1 = CombinedEigenfunction(
        triangle,
        [
            StandaloneVertexEigenfunction(triangle, 1),
            StandaloneVertexEigenfunction(triangle, 2),
            StandaloneVertexEigenfunction(triangle, 3),
            StandaloneInteriorEigenfunction(triangle),
        ],
        us_to_boundary = [
            setdiff(boundaries(triangle), (2, 3)),
            setdiff(boundaries(triangle), (1, 3)),
            setdiff(boundaries(triangle), (1, 2)),
            boundaries(triangle),
        ],
    )

    u2 = CombinedEigenfunction(
        triangle,
        [
            StandaloneVertexEigenfunction{Float64,fmpq}(triangle, 1),
            StandaloneVertexEigenfunction{Float64,fmpq}(triangle, 2),
            StandaloneVertexEigenfunction{Float64,fmpq}(triangle, 3),
            StandaloneInteriorEigenfunction(triangle),
        ],
        us_to_boundary = [
            setdiff(boundaries(triangle), (2, 3)),
            setdiff(boundaries(triangle), (1, 3)),
            setdiff(boundaries(triangle), (1, 2)),
            boundaries(triangle),
        ],
    )

    u3 = CombinedEigenfunction(
        triangle,
        [
            StandaloneLightningEigenfunction(triangle, 1),
            StandaloneLightningEigenfunction(triangle, 2),
            StandaloneLightningEigenfunction(triangle, 3),
            StandaloneInteriorEigenfunction(triangle),
        ],
    )

    u4 = CombinedEigenfunction(
        triangle,
        [
            StandaloneLightningEigenfunction{Float64,fmpq}(triangle, 1),
            StandaloneLightningEigenfunction{Float64,fmpq}(triangle, 2),
            StandaloneLightningEigenfunction{Float64,fmpq}(triangle, 3),
            StandaloneInteriorEigenfunction(triangle),
        ],
    )

    us = [u1, u2, u3, u4]
    rigorous_enclosures = [true, false, true, false]

    Ns = 8:8:24
    interval = setunion(triangle.parent.((50, 55))...)
    λs = hcat(
        [
            @time iteratemps(
                triangle,
                u,
                interval,
                Ns,
                optim_prec_final = 40,
                optim_prec_linear = true,
                show_trace = true,
                rigorous_norm = false; # Not implemented
                rigorous_enclosure,
            )[1] for (u, rigorous_enclosure) in zip(us, rigorous_enclosures)
        ]...,
    )

    for i in eachindex(Ns)
        @test overlaps(λs[i, 1], λs[i, 2])
    end

    println()
end
