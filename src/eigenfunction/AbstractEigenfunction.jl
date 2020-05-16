"""
    coefficients(u::AbstractEigenfunction)
> Return the coefficients of the eigenfunction.
"""
function coefficients(u::AbstractEigenfunction)
    u.coefficients
end

"""
    set_eigenfunction!(u::AbstractEigenfunction, coefficients)
> Set the coefficients for the expansion of the eigenfunction.
"""
function set_eigenfunction!(u::AbstractEigenfunction,
                            coefficients::Vector)
    resize!(u.coefficients, length(coefficients))
    copy!(u.coefficients, u.domain.parent.(coefficients))
end

"""
    set_domain!(u::AbstractEigenfunction, domain::AbstractDomain)
> Set the domain for the eigenfunction.
"""
function set_domain!(u::AbstractEigenfunction,
                     domain::AbstractDomain)
    u.domain = domain
    u
end

"""
    recompute!(u::AbstractEigenfunction)
> Recompute the values for the eigenfunction. This method in unsafe in
  the way that it overwrites previous values and therefore only
  handles default options.
"""
function recompute!(u::AbstractEigenfunction)
    u
end

"""
    active_boundaries(u::AbstractEigenfunction)
> Return the boundaries which are active for the current
  eigenfunction. The eigenfunction is guaranteed to be identically
  equal to zero on the inactive boundaries.
"""
function active_boundaries(u::AbstractEigenfunction)
    throw(ErrorException("active_boundaries not implemented for eigenfunction of type $(typeof(u))"))
end

"""
    u(point, λ::arb, k::Integer; notransform::Bool = false)
> Evaluate the k-th basis function for the eigenfunction at the point.

  Depending on the type of eigenfunction the type of the point varies
"""
function (u::AbstractEigenfunction)(point,
                                    λ::arb,
                                    k::Integer;
                                    boundary = nothing,
                                    notransform::Bool = false)
    throw(ErrorException("evaluation not implemented for eigenfunction of type $(typeof(u))"))
end

function norm(u::AbstractEigenfunction,
              λ::arb;
              numpoints::Int = 1000)
    @error "using a non-rigorous implementation of norm for $(typeof(u))"
    interior = interior_points(u.domain, numpoints)
    sqrt(area(u.domain)*sum(abs(u(i, λ))^2 for i in interior)/length(interior))
end
