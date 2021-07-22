@testset "example domains" begin
    @testset "decagon" begin
        println(
            "Computing eigenvalues for the decagon with different types of eigenfunctions",
        )
        parent = RealField(128)
        λ = parent(6.21200099234) # Computed by Betcke and Trefethen
        interval = setunion(λ - 1, λ + 1)

        for (lightning, linked) in
            [(false, false), (false, true), (true, false), (true, true)]
            domain, u = MethodOfParticularSolutions.example_domain_ngon(
                10,
                parent;
                lightning,
                linked,
            )

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

    @testset "Payne counterexample" begin
        # The main thing we check is that it runs correctly for the
        # different types of symmetry classes. We run it for the 1st,
        # 3rd and 4th eigenvalue. The only thing we check is that the
        # computed enclosures are finite, so not the correctness of
        # them.
        parent = RealField(64)

        N = 16
        symmetry_classes = [1, 2, 3]
        λs = [31.0432, 63.7259, 63.7259]

        for (symmetry_class, λ) in zip(symmetry_classes, λs)
            domain, u = MethodOfParticularSolutions.example_domain_goal_v1(
                27,
                11,
                6,
                parent,
                T = Float64;
                symmetry_class,
            )

            mps!(u, domain, λ, λ, N, num_boundary = 8 * N)

            enclosure = MethodOfParticularSolutions.enclose_eigenvalue_approx(
                domain,
                u,
                parent(λ),
                max_numpoints = 4 * 8 * N,
            )

            @test isfinite(enclosure)
        end
    end
end
