"""
    sigma_matrix(λ::T, domain, u, N; num_boundary, num_interior)

Return the matrix which is used for computing `σ(λ)` in the MPS.

The matrix is of size `(num_boundary + num_interior, N)`. The element
type of the matrix is `float(T)`.

The first `num_boundary` rows consists of the first `N` basis
functions of `u` evaluated on boundary points of the domain. The
remaining `num_interior` rows corresponds to the basis functions
evaluated on interior points.
"""
function sigma_matrix(
    λ::T,
    domain::AbstractDomain,
    u::AbstractEigenfunction,
    N::Integer;
    num_boundary::Integer = 2N,
    num_interior::Integer = 2N,
) where {T}
    @timeit_debug "points" begin
        # Compute boundary and interior points of the domain
        boundary, boundary_index = boundary_points(domain, u, num_boundary)
        interior = interior_points(domain, num_interior)
    end

    @timeit_debug "matrix" begin
        # Currently eigenfunctions only allow evaluation with arbs, so
        # convert λ to an arb
        λ = domain.parent(λ)

        A = Array{float(T)}(undef, num_boundary + num_interior, N)

        # TODO: Allow evaluation on ranges of indices, to be able to
        # reuse computations similar to several basis functions.
        for i in 1:num_boundary
            for k in 1:N
                A[i, k] = u(boundary[i], λ, k, boundary = boundary_index[i])
            end
        end
        for i in 1:num_interior
            for k in 1:N
                A[num_boundary + i, k] = u(interior[i], λ, k)
            end
        end
    end

    return A
end

"""
    sigma(λ, domain, u, N; num_boundary, num_interior)

Return `σ(λ)`.

Uses [`sigma_matrix`](@ref) for computing the integral which is used,
see that documentation for more details about the arguments.
"""
function sigma(
    λ::T,
    domain::AbstractDomain,
    u::AbstractEigenfunction,
    N::Integer;
    num_boundary::Integer = 2N,
    num_interior::Integer = 2N,
) where {T}
    @timeit_debug "sigma" begin
        # Compute the matrix A
        A = sigma_matrix(
            λ,
            domain,
            u,
            N,
            num_boundary = num_boundary,
            num_interior = num_interior,
        )

        # Compute a QR factorization of A
        @timeit_debug "qr" Q = Matrix(qr(A).Q)

        # Compute the smallest singular value of the top part of Q,
        # corresponding to the boundary points
        @timeit_debug "svd" value = svdvals(Q[1:num_boundary, :])[end]
    end

    return value
end

"""
    sigma_coefficients(λ, domain, u, N; num_boundary, num_interior)

Return the coefficients corresponding to `σ(λ)`.
"""
function sigma_coefficients(
    λ::T,
    domain::AbstractDomain,
    eigenfunction::AbstractEigenfunction,
    N::Integer;
    num_boundary::Integer = 2N,
    num_interior::Integer = 2N,
) where {T}
    @timeit_debug "sigma_coefficients" begin
        # Compute the matrix A
        A = sigma_matrix(
            λ,
            domain,
            eigenfunction,
            N,
            num_boundary = num_boundary,
            num_interior = num_interior,
        )

        @timeit_debug "qr" begin
            # Compute a QR factorization of A
            q = qr(A)
            Q = Matrix(q.Q)
        end

        # Compute the right singular vector corresponding to the smallest
        # singular value

        # We create a copy of Q to be able to call svd! directly,
        # otherwise it prints a warning about alg keyword being ignored.
        @timeit_debug "svd" v = svd!(
            LinearAlgebra.copy_oftype(
                Q[1:num_boundary, :],
                LinearAlgebra.eigtype(eltype(A))
            )
        ).V[:, end]

        # Compute the coefficients
        coefficients = q \ (Q*v)
    end

    return coefficients
end
