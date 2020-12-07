function IntersectedDomain{T}(
    domain::IntersectedDomain{T},
    parent::ArbField,
) where {T<:AbstractPlanarDomain}
    return IntersectedDomain(
        T(domain.exterior, parent),
        [typeof(d)(d, parent) for d in domain.interiors],
    )
end

IntersectedDomain(exterior::AbstractPlanarDomain, interior::AbstractPlanarDomain) =
    IntersectedDomain(exterior, [interior])

function Base.show(io::IO, domain::IntersectedDomain)
    println(io, "Intersected domain")
    println(io, "Exterior: $(domain.exterior)")
    println(io, "Interiors:")
    for d in domain.interiors
        println(d)
    end
end

function Base.getproperty(domain::IntersectedDomain, name::Symbol)
    if name == :parent
        return domain.exterior.parent
    else
        return getfield(domain, name)
    end
end


exterior_boundaries(domain::IntersectedDomain) = boundaries(domain.exterior)

function interior_boundaries(domain::IntersectedDomain)
    boundaries_exterior = boundaries(domain.exterior)
    boundaries_interior = boundaries.(domain.interiors)
    return (1:sum(length.(boundaries_interior))) .+ maximum(boundaries_exterior)
end

boundaries(domain::IntersectedDomain) =
    1:(length(exterior_boundaries(domain)) + length(interior_boundaries(domain)))

function get_domain_and_boundary(domain::IntersectedDomain, i::Integer)
    boundaries_exterior = boundaries(domain.exterior)
    if i ∈ boundaries_exterior
        return domain.exterior, i
    end

    i_local = i - maximum(boundaries_exterior)
    for d in domain.interiors
        bs = boundaries(d)
        if i_local ∈ bs
            return (d, i_local)
        end
        i_local -= maximum(bs)
    end

    throw(ArgumentError("attempt to get boundary $i from a $(typeof(domain))"))
end

function vertex(domain::IntersectedDomain, i::Integer)
    d, j = get_domain_and_boundary(domain, i)
    return vertex(d, j)
end

vertices(domain::IntersectedDomain) = [
    vertices(domain.exterior)...,
    [vertices.(domain.interiors)...]...,
]

center(domain::IntersectedDomain) = center(domain.exterior)

area(domain::IntersectedDomain) = area(domain.exterior) - sum(area.(domain.interiors))

Base.in(xy, domain::IntersectedDomain) =
    xy ∈ domain.exterior && !any(d -> xy ∈ d, domain.interiors)

boundary_parameterization(t, domain::IntersectedDomain, i::Integer) =
    boundary_parameterization(t, get_domain_and_boundary(domain, i)...)

boundary_points(domain::IntersectedDomain, i::Integer, n::Integer) =
    boundary_points(get_domain_and_boundary(domain, i)..., n)

# TODO: Exclude points landing in the interior domain
function interior_points(domain::IntersectedDomain, n::Integer; rng = MersenneTwister(42))
    # We generate points in the interior of domain.exterior and remove
    # those which are also in domain.interior. The issue is that we
    # then don't know how many interior points we need to generate.
    gotten = 0
    points = SVector{2, arb}[]
    while gotten < n
        new_points = Iterators.take(
            filter(
                xy -> !any(d -> xy ∈ d, domain.interiors),
                interior_points(domain.exterior, 2(n - gotten), rng = rng),
            ),
            n - gotten
        )

        gotten += length(new_points)
        append!(points, new_points)
    end

    return points
end
