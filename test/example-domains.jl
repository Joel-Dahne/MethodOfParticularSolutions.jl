@testset "example domains" begin
    println("Computing eigenvalues for the decagon with different types of eigenfunctions")
    parent = RealField(128)
    λ = parent(6.21200099234) # Computed by Betcke and Trefethen
    interval = setunion(λ - 1, λ + 1)

    for (lightning, linked) in [(false, false), (false, true), (true, false), (true, true)]
        domain, u =
            MethodOfParticularSolutions.example_domain_ngon(10, parent; lightning, linked)

        λs, us = iteratemps(
            domain,
            u,
            interval,
            ifelse(linked, 2:2:4, 8:8:16),
            num_boundary_factor = ifelse(lightning && linked, 8, 2),
            optim_prec_final = 20,
            show_trace = true,
            rigorous_enclosure = true,
            rigorous_norm = false, # Not implemented
        )

        @test all(contains(λ), λs)
    end

    println()
end
