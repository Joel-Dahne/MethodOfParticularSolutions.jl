Polygon{S,T}(
    domain::Polygon{S,T};
    parent::Union{ArbField,Nothing} = domain.parent,
) where {S,T} = Polygon{S,T}(copy(domain.angles), copy(domain.vertices); parent)

function Base.show(io::IO, domain::Polygon)
    print(io, "Polygon with $(length(boundaries(domain))) sides")
    if !isnothing(domain.parent)
        print(io, " and $(precision(domain.parent)) bits of precision")
    end
end

vertexindices(domain::Polygon) = 1:length(domain.angles)
boundaries(domain::Polygon) = 1:length(domain.angles)

angle_raw(domain::Polygon, i::Integer) = domain.angles[i]

vertex(domain::Polygon, i::Integer) = domain.vertices[i]
vertices(domain::Polygon) = domain.vertices

function orientation_raw(
    domain::Polygon{S,<:AbstractFloat},
    i::Integer;
    reversed = false,
) where {S}
    # Compute the orientation for the first vertex and then use the
    # angles to get the remaining orientations. Currently we could
    # compute the orientation directly, but this will be better if we
    # switch to storing an exact orientation.
    # FIXME: This is not accurate if v[2] contains zero but is not
    # exactly zero
    v = vertex(domain, 2) - vertex(domain, 1)
    res = atan(v[2], v[1])

    res += (i - 1) * oftype(res, π)
    for j = 2:i
        res -= angle_raw(domain, j)
    end

    if reversed
        res = 2oftype(res, π) - angle_raw(domain, i) - res
    end

    return res
end

function orientation_raw(domain::Polygon{arb,arb}, i::Integer; reversed = false)
    # Compute the orientation for the first vertex and then use the
    # angles to get the remaining orientations. Currently we could
    # compute the orientation directly, but this will be better if we
    # switch to storing an exact orientation.
    # FIXME: This will not be accurate if v[2] contains zero but is
    # not exactly zero
    v = vertex(domain, 2) - vertex(domain, 1)
    res = atan(v[2], v[1])

    res += (i - 1) * domain.parent(π)
    for j = 2:i
        res -= angle_raw(domain, j)
    end

    if reversed
        res = 2domain.parent(π) - angle_raw(domain, i) - res
    end

    return res
end

function orientation(domain::Polygon{S,<:Rational}, i::Integer; reversed = false) where {S}
    # Compute the orientation for the first vertex and then use the
    # angles to get the remaining orientations. Currently we could
    # compute the orientation directly, but this will be better if we
    # switch to storing an exact orientation.
    # FIXME: This is only needed because we can't implement an
    # orientation_raw method for rational angles
    v = vertex(domain, 2) - vertex(domain, 1)
    res = atan(v[2], v[1])

    res += (i - 1) * oftype(res, π)
    for j = 2:i
        res -= angle(domain, j)
    end

    if reversed
        res = 2oftype(res, π) - angle(domain, i) - res
    end

    return res
end

function orientation(domain::Polygon{arb,fmpq}, i::Integer; reversed = false)
    # Compute the orientation for the first vertex and then use the
    # angles to get the remaining orientations. Currently we could
    # compute the orientation directly, but this will be better if we
    # switch to storing an exact orientation.
    # FIXME: This will not be accurate if v[2] contains zero but is
    # not exactly zero.
    # FIXME: This is only needed because we can't implement an
    # orientation_raw method for rational angles
    v = vertex(domain, 2) - vertex(domain, 1)
    res = atan(v[2], v[1])

    res += (i - 1) * domain.parent(π)
    for j = 2:i
        res -= angle(domain, j)
    end

    if reversed
        res = 2domain.parent(π) - angle(domain, i) - res
    end

    return res
end

function area(domain::Polygon{S,T}) where {S,T}
    # https://en.wikipedia.org/wiki/Polygon#Area
    A = zero(S)
    vs = vertices(domain)
    for i in eachindex(vs)
        A +=
            vs[i][1] * vs[mod1(i + 1, length(vs))][2] -
            vs[mod1(i + 1, length(vs))][1] * vs[i][2]
    end

    return abs(A) / 2
end

function area(domain::Polygon{arb,T}) where {T}
    # https://en.wikipedia.org/wiki/Polygon#Area
    A = domain.parent(0)
    vs = vertices(domain)
    for i in eachindex(vs)
        A +=
            vs[i][1] * vs[mod1(i + 1, length(vs))][2] -
            vs[mod1(i + 1, length(vs))][1] * vs[i][2]
    end

    return abs(A) / 2
end

center(domain::Polygon) = sum(vertices(domain)) / length(boundaries(domain))

function Base.in(xy, domain::Polygon)
    vs = vertices(domain)

    isleft(p0, p1, p2) =
        (p1[1] - p0[1]) * (p2[2] - p0[2]) - (p2[1] - p0[1]) * (p1[2] - p0[2])

    wn = 0
    for i in eachindex(vs)
        vi = vs[i]
        vip1 = vs[mod1(i + 1, length(vs))]
        if vi[2] <= xy[2]
            if !(vip1[2] <= xy[2])
                if isleft(vi, vip1, xy) > 0
                    wn += 1
                end
            end
        else
            if vip1[2] <= xy[2]
                if isleft(vi, vip1, xy) < 0
                    wn -= 1
                end
            end
        end
    end

    return !iszero(wn)
end

function boundary_parameterization(t, domain::Polygon, i::Integer)
    i ∈ boundaries(domain) ||
        throw(ArgumentError("attempt to get vertex $i from a $(typeof(domain))"))
    v = vertex(domain, mod1(i, length(boundaries(domain))))
    w = vertex(domain, mod1(i + 1, length(boundaries(domain))))

    return v .+ t .* (w - v)
end

function boundary_points(
    domain::Polygon{S},
    i::Integer,
    n::Integer;
    distribution = :chebyshev,
) where {S}
    i ∈ boundaries(domain) ||
        throw(ArgumentError("attempt to get vertex $i from a $(typeof(domain))"))

    v = vertex(domain, mod1(i, length(boundaries(domain))))
    w = vertex(domain, mod1(i + 1, length(boundaries(domain))))

    if distribution == :linear
        if S == arb
            f = j -> domain.parent(j // (n + 1))
        else
            f = j -> convert(S, j // (n + 1))
        end
    elseif distribution == :chebyshev
        if S == arb
            f = j -> (cospi(domain.parent((2j - 1) // 2n)) + 1) / 2
        else
            f = j -> (cospi(convert(S, (2j - 1) // 2n)) + 1) / 2
        end
    elseif distribution == :exponential
        m = (n + 1) // 2
        if S == arb
            f = j -> begin
                d = exp(domain.parent(m - j))
                if j < m
                    return d / 2
                elseif j > m
                    return 1 - inv(d) / 2
                else
                    return domain.parent(0.5)
                end
            end
        else
            f = j -> begin
                d = exp(-convert(S, m - j))
                if j < m
                    return d / 2
                elseif j > m
                    return 1 - 1 / 2d
                else
                    return convert(S, 1 // 2)
                end
            end
        end
    elseif distribution == :root_exponential
        m = (n + 1) // 2
        if S == arb
            f =
                j -> begin
                    d = exp(-4(sqrt(domain.parent(m)) - sqrt(domain.parent(min(j, n - j + 1)))))
                    if j < m
                        return d / 2
                    elseif j > m
                        return 1 - d / 2
                    else
                        return domain.parent(0.5)
                    end
                end
        else
            f =
                j -> begin
                    d = exp(-4(sqrt(convert(S, m)) - sqrt(convert(S, min(j, n - j + 1)))))
                    if j < m
                        return d / 2
                    elseif j > m
                        return 1 - d / 2
                    else
                        return convert(S, 1 // 2)
                    end
                end
        end
    else
        throw(ArgumentError("no distribution named $distribution"))
    end

    points = Vector{typeof(v)}(undef, n)
    for j = 1:n
        t = f(j)
        points[j] = v .+ t .* (w - v)
    end

    return points, fill(i, n)
end

function interior_points(
    domain::Polygon{S},
    n::Integer;
    rng = MersenneTwister(42),
) where {S}
    points = Vector{SVector{2,S}}(undef, n)
    xmin, xmax = extrema(getindex.(vertices(domain), 1))
    ymin, ymax = extrema(getindex.(vertices(domain), 2))
    for i = 1:n
        pt = SVector(xmin + rand(rng) * (xmax - xmin), ymin + rand(rng) * (ymax - ymin))
        while !(pt ∈ domain)
            pt = SVector(xmin + rand(rng) * (xmax - xmin), ymin + rand(rng) * (ymax - ymin))
        end

        points[i] = pt
    end

    return points
end
