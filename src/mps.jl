function mps(domain::AbstractDomain,
             eigenfunction::AbstractEigenfunction,
             enclosure::arb,
             N::Integer = 8;
             num_boundary = 2N,
             num_interior = 2N,
             store_trace = false,
             show_trace = false,
             extended_trace = false,
             optim_rel_tol = eps(domain.parent),
             optim_abs_tol = eps(domain.parent),
             optim_store_trace = false,
             optim_show_trace = false,
             optim_extended_trace = false,
             enclose_store_trace = false,
             enclose_show_trace = false,
             enclose_extended_trace = false)
    # Setup

    # Compute minimum of σ(λ)
    σ = λ -> sigma(λ, domain, eigenfunction, N,
                   num_boundary = num_boundary,
                   num_interior = num_interior)
    res = optimize(σ,
                   getinterval(BigFloat, enclosure)...,
                   rel_tol = BigFloat(optim_rel_tol),
                   abs_tol = BigFloat(optim_abs_tol),
                   store_trace = optim_store_trace,
                   show_trace = optim_show_trace,
                   extended_trace = optim_extended_trace)

    if !Optim.converged(res)
        @warn "Failed to compute minimum of σ(λ)"
    end

    λ = res.minimizer

    # Compute the eigenfunction corresponding to the computed
    # minimum of σ(λ).
    coefficients = sigma_coefficients(λ, domain, eigenfunction, N,
                                      num_boundary = num_boundary,
                                      num_interior = num_interior)
    set_eigenfunction!(eigenfunction, coefficients)

    # Compute the enclosure of the eigenvalue
    λ = enclose_eigenvalue(domain,
                           eigenfunction,
                           domain.parent(λ),
                           store_trace = enclose_store_trace,
                           show_trace = enclose_show_trace,
                           extended_trace = enclose_extended_trace)

    return λ, eigenfunction
end
