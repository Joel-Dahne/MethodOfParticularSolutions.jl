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
    boundary_points(
        domain::AbstractDomain,
        u::AbstractEigenfunction,
        N::Integer,
        n::Integer;
        distribution,
)

Return `n` points from boundaries of the domain on which the
eigenfunction is not identically zero, together with information about
which boundary they are from.

The boundaries where the egienfunction is not identically equal to
zero on is determined by `active_boundaries(domain, u)`.

`N` should be set to the number of expansion terms used for `u`. It's
not used in the default implementation but some type of eigenfunctions
need access to this to compute good boundary points.
"""
function boundary_points(
    domain::AbstractDomain,
    u::AbstractEigenfunction,
    ::Integer,
    n::Integer;
    distribution = ifelse(domain isa AbstractPlanarDomain, :chebyshev, :linear),
)
    active = active_boundaries(domain, u)
    m = length(active)
    res = [
        boundary_points(
            domain,
            active[i],
            div(n, m) + ifelse(n % m >= i, 1, 0);
            distribution,
        ) for i in eachindex(active)
    ]
    return vcat(getindex.(res, 1)...), vcat(getindex.(res, 2)...)
end

function boundary_points(
    domain::AbstractDomain,
    u::CombinedEigenfunction,
    ::Integer,
    n::Integer;
    distribution = ifelse(domain isa AbstractPlanarDomain, :root_exponential, :linear),
)
    active = active_boundaries(domain, u)
    m = length(active)
    per_boundary = fill(div(n, m), m) + [n % m >= i for i = 1:m]
    res = [
        getindex.(
            boundary_points(
                domain,
                active[i],
                ifelse(
                    active[i] ∈ u.even_boundaries,
                    2per_boundary[i] - 1,
                    per_boundary[i],
                );
                distribution,
            ),
            Ref(1:per_boundary[i]),
        ) for i in eachindex(active)
    ]
    return vcat(getindex.(res, 1)...), vcat(getindex.(res, 2)...)
end

"""
    interior_points(domain::AbstractDomain, n::Integer; rng = MersenneTwister(42))

Return n points taken uniformly at random from the interior of the
domain.
"""
interior_points(domain::AbstractDomain, n::Integer; rng = MersenneTwister(42))
