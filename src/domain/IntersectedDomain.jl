function Base.show(io::IO, domain::IntersectedDomain)
    println(io, "Intersected domain")
    print(io, "Exterior: $(domain.exterior)")
    print(io, "Interior: $(domain.interior)")
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
    boundaries_interior = boundaries(domain.interior)
    return boundaries_interior .+ maximum(boundaries_exterior)
end

function boundaries(domain::IntersectedDomain)
    boundaries_exterior = boundaries(domain.exterior)
    boundaries_interior = boundaries(domain.interior)
    return union(boundaries_exterior, boundaries_interior .+ maximum(boundaries_exterior))
end

function get_domain_and_boundary(domain::IntersectedDomain, i::Integer)
    boundaries_exterior = boundaries(domain.exterior)
    if i ∈ boundaries_exterior
        return domain.exterior, i
    end

    boundaries_interior = boundaries(domain.interior)
    if i - maximum(boundaries_exterior) ∈ boundaries_interior
        return domain.interior, i - maximum(boundaries_exterior)
    end

    throw(ArgumentError("attempt to get vertex $i from a $(typeof(domain))"))
end

center(domain::IntersectedDomain) = center(domain.exterior)

area(domain::IntersectedDomain) = area(domain.exterior) - area(domain.interior)

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
                xy -> !(xy ∈ domain.interior),
                interior_points(domain.exterior, 2(n - gotten), rng = rng),
            ),
            n - gotten
        )

        gotten += length(new_points)
        append!(points, new_points)
    end

    return points
end
