TransformedDomain{S,T,U}(
    domain::TransformedDomain{S,T,U};
    parent::Union{ArbField,Nothing} = domain.parent,
) where {S,T,U} = TransformedDomain(
    S(domain.original; parent),
    domain.rotation,
    domain.scaling,
    domain.translation,
)

function Base.show(io::IO, domain::TransformedDomain{T}) where {T}
    println(io, "Transformed Domain")
    if T == arb
        println(io, "Rotation: $(domain.rotation)")
    else
        println(
            io,
            "Rotation: $(numerator(domain.rotation))π/$(denominator(domain.rotation))",
        )
    end
    println(io, "Scaling: $(domain.scaling)")
    println(io, "Translation: $(domain.translation)")
    print(io, "Domain: $(domain.original)")
end

function Base.getproperty(domain::TransformedDomain{S,T}, name::Symbol) where {S,T}
    if name == :parent
        return domain.original.parent
    elseif name == :map
        if T == fmpq
            s, c = sincospi(domain.rotation, domain.parent)
        elseif has_rational_angles(domain)
            s, c = sincospi(convert(S, domain.rotation))
        else
            s, c = sincos(domain.rotation)
        end
        M = SMatrix{2,2}(c, s, -s, c)
        return xy -> domain.scaling .* M * xy + domain.translation
    elseif name == :invmap
        if T == fmpq
            s, c = sincospi(-domain.rotation, domain.parent)
        elseif has_rational_angles(domain)
            s, c = sincospi(convert(S, -domain.rotation))
        else
            s, c = sincos(-domain.rotation)
        end
        M = SMatrix{2,2}(c, s, -s, c)
        return xy -> M * (xy - domain.translation) ./ domain.scaling
    else
        return getfield(domain, name)
    end
end

vertexindices(domain::TransformedDomain) = vertexindices(domain.original)
boundaries(domain::TransformedDomain) = boundaries(domain.original)

angle_raw(domain::TransformedDomain, i::Integer) = angle_raw(domain.original, i)

vertex(domain::TransformedDomain, i::Integer) = domain.map(vertex(domain.original, i))

orientation_raw(domain::TransformedDomain, i::Integer; reversed = false) =
    orientation_raw(domain.original, i; reversed) +
    ifelse(reversed, -domain.rotation, domain.rotation)

center(domain::TransformedDomain) = domain.map(center(domain.original))

area(domain::TransformedDomain) = domain.scaling * area(domain.original)

Base.in(xy, domain::TransformedDomain) = domain.invmap(xy) ∈ domain.original

boundary_parameterization(t, domain::TransformedDomain, i::Integer) =
    domain.map(boundary_parameterization(t, domain.original, i))

function boundary_points(
    domain::TransformedDomain,
    i::Integer,
    n::Integer;
    distribution = :chebyshev,
)
    pts, idx = boundary_points(domain.original, i, n; distribution)
    return domain.map.(pts), idx
end

interior_points(domain::TransformedDomain, n::Integer; rng = MersenneTwister(42)) =
    domain.map.(interior_points(domain.original, n, rng = rng))
