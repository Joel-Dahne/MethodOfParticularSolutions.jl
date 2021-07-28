"""
    has_rational_angles(u::AbstractEigenfunction)

Return `true` if the angles of `u` are represented as rational
multiples of `π`, return `false` otherwise.

Notice that this only considers how the angles are represented, not
the actual values. The domain could have angles which are rational
multiples of `π` but are stored as floating points.
"""
has_rational_angles(::AbstractEigenfunction{T,<:Union{Rational,fmpq}}) where {T} = true
has_rational_angles(::AbstractEigenfunction{T,<:Union{AbstractFloat,arb}}) where {T} = false

"""
    coefficients(u::AbstractEigenfunction)

Return the coefficients of the eigenfunction.
"""
coefficients(u::AbstractEigenfunction) = u.coefficients

"""
    set_eigenfunction!(u::AbstractEigenfunction, coefficients)

Set the coefficients for the expansion of the eigenfunction.
"""
function set_eigenfunction!(u::AbstractEigenfunction{arb}, coefficients::Vector)
    resize!(u.coefficients, length(coefficients))
    copy!(u.coefficients, u.parent.(coefficients))
    return u
end

function set_eigenfunction!(u::AbstractEigenfunction, coefficients::Vector)
    resize!(u.coefficients, length(coefficients))
    copy!(u.coefficients, coefficients)
    return u
end

"""
    set_domain!(u::AbstractEigenfunction, domain::AbstractDomain)

If the eigenfunction stores the domain update it to the given domain.
"""
function set_domain!(u::AbstractEigenfunction, domain::AbstractDomain)
    if hasfield(typeof(u), :domain)
        u.domain = domain
    end
    return u
end

"""
    recompute!(u::AbstractEigenfunction)

Recompute the values for the eigenfunction.

Some type of eigenfunctions have precomputed values which might need
to be recomputed if the precision is increased, then this method can
be used. However this is mostly older versions of eigenfunctions. This
method in unsafe in the way that it overwrites previous values and
therefore only handles default options.
"""
recompute!(u::AbstractEigenfunction) = u

"""
    active_boundaries(domain::AbstractDomain, u::AbstractEigenfunction)

Return the boundaries on `domain` which are active for `u`. The
eigenfunction is guaranteed to be identically equal to zero on the
inactive boundaries.

Some eigenfunctions are constructed for a specific domain, given by
`u.doman`. In this case it checks if `domain === u.domain` (i.e. if
they are the exact same object in Julia), if that is true then it uses
information `u` has about its domain to determine which boundaries are
active. If `domain !== u.domain` the it just returns all boundaries of
`domain`.

Many eigenfunctions are not tied to a domain directly, in that case it
will in general just return all boundaries of `domain`.

One of the reason it works like this is that it's hard to determine if
two domains are exactly the same since we generally work with ball
arithmetic. It's also a workaround to have the old version which only
relied on `u.domain` and the new version where `u` might not have a
specified domain work together.
"""
active_boundaries(domain::AbstractDomain, u::AbstractEigenfunction) = boundaries(domain)

"""
    u(point, λ, k::Integer; boundary = nothing, notransform::Bool = false)

Evaluate the `k`-th basis function for the eigenfunction at the point.

Depending on the type of eigenfunction the type of the point varies
"""
function (u::AbstractEigenfunction)(
    point,
    λ,
    k::Integer;
    boundary = nothing,
    notransform::Bool = false,
) end

function norm(
    domain::AbstractDomain,
    u::AbstractEigenfunction,
    λ::arb;
    numpoints::Int = 1000,
    warn::Bool = true,
)
    warn && @warn "using a non-rigorous implementation of norm for $(typeof(u))"
    interior = interior_points(domain, numpoints)
    res = similar(interior, arb)
    Threads.@threads for i in eachindex(interior)
        res[i] = abs(u(interior[i], λ))^2
    end
    return sqrt(area(domain) * sum(res) / length(interior))
end
