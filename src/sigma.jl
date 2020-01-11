"""
    sigma(λ, ...)
> Return σ(λ).
"""
function sigma(λ,
               domain::AbstractDomain,
               eigenfunction::AbstractEigenfunction,
               N::Integer;
               num_boundary::Integer = 2N,
               num_interior::Integer = 2N)
    # Compute boundary and interior points of the domain
    boundary = boundary_points(domain, eigenfunction, num_boundary)
    interior = interior_points(domain, num_interior)

    # Evaluate the basis of the eigenfunction on the points
    A_B = BigFloat.(eigenfunction.(boundary, domain.parent(λ), (1:N)'))
    A_I = BigFloat.(eigenfunction.(interior, domain.parent(λ), (1:N)'))

    A = [A_B; A_I]

    # Compute a QR factorization of A
    Q = Matrix(qr(A).Q)

    # Return the smallest singular value of the top part of Q
    # corresponding to the interior points
    svdvals(Q[1:num_boundary, :])[end]
end

"""
    sigma_coefficients(λ, )
> Return the coefficients corresponding to σ(λ).
"""
function sigma_coefficients(λ,
                            domain::AbstractDomain,
                            eigenfunction::AbstractEigenfunction,
                            N::Integer;
                            num_boundary::Integer = 2N,
                            num_interior::Integer = 2N)
    # Compute boundary and interior points of the domain
    boundary = boundary_points(domain, eigenfunction, num_boundary)
    interior = interior_points(domain, num_interior)

    # Evaluate the basis of the eigenfunction on the points
    A_B = BigFloat.(eigenfunction.(boundary, domain.parent(λ), (1:N)'))
    A_I = BigFloat.(eigenfunction.(interior, domain.parent(λ), (1:N)'))

    A = [A_B; A_I]

    # Compute a QR factorization of A
    q = qr(A)
    Q = Matrix(q.Q)

    # Compute the right singular vector corresponding to the smallest
    # singular value
    v = svd(Q[1:num_boundary, :]).V[:, end]

    # Compute and return the coefficients
    coefficients = q \ (Q*v)
end
