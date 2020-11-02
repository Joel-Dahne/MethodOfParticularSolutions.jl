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
    active_boundaries(domain::AbstractDomain, u::AbstractEigenfunction)

Return the boundaries on `domain` which are active for `u`. The
eigenfunction is guaranteed to be identically equal to zero on the
inactive boundaries.

Most eigenfunctions are constructed for a specific domain, given by
`u.doman`. In this case it checks if `domain === u.domain` (i.e. if
they are the exact same object in Julia), if that is true then it uses
information `u` has about it's domain to determine which boundaries
are active. If `domain !== u.domain` the it just returns all
boundaries of `domain`.

Some eigenfunctions are not tied to a domain directly, in that case it
will in general just return all boundaries of `domain`.

One of the reason it works like this is that it's hard to determine of
two domains are exactly the same since we generally work with ball
arithmetic. It's also a workaround to have the old version which only
relied on `u.domain` and the new version where `u` might not have a
specified domain work together.
"""
active_boundaries(domain::AbstractDomain, u::AbstractEigenfunction) = boundaries(domain)

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

function norm(domain::AbstractDomain,
              u::AbstractEigenfunction,
              λ::arb;
              numpoints::Int = 1000)
    @error "using a non-rigorous implementation of norm for $(typeof(u))"
    interior = interior_points(domain, numpoints)
    sqrt(area(domain)*sum(abs(u(i, λ))^2 for i in interior)/length(interior))
end
