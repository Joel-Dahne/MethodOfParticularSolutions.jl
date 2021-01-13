"""
    boundaries(domain::AbstractDomain)

Return an index set giving the boundaries of the domain.
"""
boundaries(domain::AbstractDomain)

"""
    area(domain::AbstractDomain)

Return the area of the domain.
"""
area(domain::AbstractDomain)

"""
    boundary_parameterization(t, domain::AbstractDomain, i::Integer)

Compute a parameterization of boundary number `i`.

The parameterization goes from t = 0 to t = 1. The enumeration of
the boundaries depends on the type of domain used.
"""
boundary_parameterization(t, domain::AbstractDomain, i::Integer)

"""
    boundary_points(domain::AbstractDomain, i::Integer, n::Integer; distribution)

Return `n` points from boundary number `i`, together with information
about which boundary they are from.

The enumeration of the boundary depends on the type of domain used.
The distribution of the points depends on `distribution`, the default
is `:chebyshev` for planar domains and `:linear` for spherical ones.
possible options depends on the domain.
"""
boundary_points(domain::AbstractDomain, i::Integer, n::Integer; distribution = :chebyshev)

"""
    boundary_points(domain::AbstractDomain, n::Integer)

Return `n` points taken from all boundaries of the domain, together
with information about which boundary they are from.
"""
function boundary_points(
    domain::AbstractDomain,
    n::Integer;
    distribution = ifelse(domain isa AbstractPlanarDomain, :chebyshev, :linear),
)
    m = length(boundaries(domain))
    res = [
        boundary_points(domain, i, div(n, m) + ifelse(n % m >= i, 1, 0); distribution)
        for i in boundaries(domain)
    ]
    return vcat(getindex.(res, 1)...), vcat(getindex.(res, 2)...)
end

"""
    boundary_points(domain::AbstractDomain, u::AbstractEigenfunction, n::Integer; distribution)

Return n points from boundaries of the domain on which the
eigenfunction is not identically zero, together with information about
which boundary they are from.

The boundaries where the egienfunction is not identically equal to
zero on is determined by `active_boundaries(domain, u)`.
"""
function boundary_points(domain::AbstractDomain,
                         u::AbstractEigenfunction,
                         n::Integer;
                         distribution = ifelse(
                             domain isa AbstractPlanarDomain,
                             :chebyshev,
                             :linear,
                         ),
                         )
    active = active_boundaries(domain, u)
    m = length(active)
    res = [
        boundary_points(domain, active[i], div(n, m) + ifelse(n % m >= i, 1, 0); distribution)
        for i in eachindex(active)
    ]
    return vcat(getindex.(res, 1)...), vcat(getindex.(res, 2)...)
end

"""
    interior_points(domain::AbstractDomain, n::Integer; rng = MersenneTwister(42))

Return n points taken uniformly at random from the interior of the
domain.
"""
interior_points(domain::AbstractDomain, n::Integer; rng = MersenneTwister(42))
