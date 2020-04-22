function sigma_matrix(λ,
                      domain::AbstractDomain,
                      eigenfunction::AbstractEigenfunction,
                      N::Integer;
                      num_boundary::Integer = 2N,
                      num_interior::Integer = 2N)
    @timeit_debug "points" begin
        # Compute boundary and interior points of the domain
        boundary, boundary_index = boundary_points(domain, eigenfunction, num_boundary)
        interior = interior_points(domain, num_interior)

    end

    @timeit_debug "matrix" begin
        A = Array{BigFloat}(undef, num_interior + num_boundary, N)
        for k in 1:N
            for i in 1:num_boundary
                A[i, k] = BigFloat(eigenfunction(boundary[i], domain.parent(λ), k, boundary = boundary_index[i]))
            end

            for i in 1:num_interior
                A[num_boundary + i, k] = BigFloat(eigenfunction(interior[i], domain.parent(λ), k))
            end
        end
    end

    A
end

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
    @timeit_debug "sigma" begin
        # Compute the matrix A
        A = sigma_matrix(λ, domain, eigenfunction, N,
                         num_boundary = num_boundary,
                         num_interior = num_interior)

        # Compute a QR factorization of A
        @timeit_debug "qr" Q = Matrix(qr(A).Q)

        # Return the smallest singular value of the top part of Q
        # corresponding to the interior points
        @timeit_debug "svd" svdvals(Q[1:num_boundary, :])[end]
    end
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
    @timeit_debug "sigma_coefficients" begin
        # Compute the matrix A
        A = sigma_matrix(λ, domain, eigenfunction, N,
                         num_boundary = num_boundary,
                         num_interior = num_interior)

        @timeit_debug "qr" begin
            # Compute a QR factorization of A
            q = qr(A)
            Q = Matrix(q.Q)
        end

        # Compute the right singular vector corresponding to the smallest
        # singular value

        # We create a copy of Q to be able to call svd! directly,
        # otherwise it prints a warning about alg keyword being ignored.
        @timeit_debug "svd" v = svd!(LinearAlgebra.copy_oftype(Q[1:num_boundary, :],
                                                               LinearAlgebra.eigtype(BigFloat))).V[:, end]

        # Compute and return the coefficients
        coefficients = q \ (Q*v)
    end
end
