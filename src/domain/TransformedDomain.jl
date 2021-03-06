function TransformedDomain{T,S}(
    domain::TransformedDomain{T,S},
    parent::ArbField,
) where {T,S}
    # TODO: Might need to recompute values here
    return TransformedDomain(
        S(domain.original, parent),
        domain.rotation,
        domain.scaling,
        domain.translation,
    )
end

function TransformedDomain(domain, rotation::T, scaling, translation) where {T}
    if T == fmpq || T <: Rational
        rotation = fmpq(rotation)
    else
        rotation = domain.parent(rotation)
    end
    scaling = domain.parent(scaling)
    translation = SVector{2,arb}(domain.parent.(translation))
    return TransformedDomain(domain, rotation, scaling, translation)
end

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

function Base.getproperty(domain::TransformedDomain{T}, name::Symbol) where {T}
    if name == :parent
        return domain.original.parent
    elseif name == :map
        if T == arb
            s, c = sincos(domain.rotation)
        elseif T == fmpq
            s, c = sincospi(domain.rotation, domain.parent)
        end
        M = SMatrix{2,2}(c, s, -s, c)
        return xy -> domain.scaling .* M * xy + domain.translation
    elseif name == :invmap
        if T == arb
            s, c = sincos(-domain.rotation)
        elseif T == fmpq
            s, c = sincospi(-domain.rotation, domain.parent)
        end
        M = SMatrix{2,2}(c, s, -s, c)
        return xy -> M * (xy - domain.translation) ./ domain.scaling
    else
        return getfield(domain, name)
    end
end

boundaries(domain::TransformedDomain) = boundaries(domain.original)

angle(domain::TransformedDomain, i::Integer) = angle(domain.original, i)
angledivπ(domain::TransformedDomain, i::Integer) = angledivπ(domain.original, i)
angles(domain::TransformedDomain) = angles(domain.original)
angledivπ(domain::TransformedDomain) = angledivπ(domain.original)

vertex(domain::TransformedDomain, i::Integer) = domain.map(vertex(domain.original, i))
vertices(domain::TransformedDomain) = domain.map.(vertices(domain.original))

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
