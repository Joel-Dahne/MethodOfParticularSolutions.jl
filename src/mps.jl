function mps(domain::AbstractDomain,
             eigenfunction::AbstractEigenfunction,
             enclosure::arb,
             N::Integer = 8)
    # Setup

    # Compute minimum of σ(λ)
    σ = λ -> sigma(λ, domain, eigenfunction, N)
    res = optimize(σ, getinterval(BigFloat, enclosure)...)

    if !Optim.converged(res)
        @warn "Failed to compute minimum of σ(λ)"
    end

    # Compute the eigenfunction corresponding to the computed
    # minimum of σ(λ).
    coefficients = sigma_coefficients(λ, domain, eigenfunction, N)
    set_eigenfunction!(eigenfunction, coefficients)

    # Return the enclosure of the eigenvalue
    enclose_eigenvalue(domain, eigenfunction)
end
