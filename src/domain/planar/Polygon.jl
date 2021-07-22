function Polygon{T}(
    domain::Polygon{T},
    parent::ArbField = domain.parent,
) where {T<:Union{arb,fmpq}}
    return Polygon(copy(domain.angles), copy(domain.vertices), parent)
end


function Polygon(
    angles::AbstractVector{T},
    vertices::AbstractVector,
    parent::ArbField,
) where {T}
    @assert length(angles) == length(vertices)
    if T == fmpq || T <: Rational
        S = fmpq
        angles = fmpq.(angles)
    else
        S = arb
        angles = parent.(angles)
    end
    vertices = [SVector{2,arb}(parent.(vertex)) for vertex in vertices]
    return Polygon{S}(angles, vertices, parent)
end


function Base.show(io::IO, domain::Polygon)
    print(
        io,
        "Polygon with $(length(boundaries(domain))) sides and $(domain.parent.prec) bits of precision",
    )
end

boundaries(domain::Polygon) = 1:length(domain.angles)

angle(domain::Polygon{fmpq}, i::Integer) = domain.parent(π) * domain.angles[i]
angle(domain::Polygon{arb}, i::Integer) = domain.angles[i]

angledivπ(domain::Polygon{fmpq}, i::Integer) = domain.angles[i]
angledivπ(domain::Polygon{arb}, i::Integer) = domain.angles[i] / domain.parent(π)

angles(domain::Polygon) = [angle(domain, i) for i in boundaries(domain)]
anglesdivπ(domain::Polygon) = [angledivπ(domain, i) for i in boundaries(domain)]

vertex(domain::Polygon, i::Integer) = domain.vertices[i]
vertices(domain::Polygon) = domain.vertices

function orientation(domain::Polygon, i::Integer; reversed = false)
    # Compute the orientation for the first vertex and then use the
    # angles to get the remaining orientations. Currently we could
    # compute the orientation directly, but this will be better if we
    # switch to storing an exact orientation.
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

function area(domain::Polygon)
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

# TODO: Add this
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

function boundary_points(domain::Polygon, i::Integer, n::Integer; distribution = :chebyshev)
    i ∈ boundaries(domain) ||
        throw(ArgumentError("attempt to get vertex $i from a $(typeof(domain))"))

    v = vertex(domain, mod1(i, length(boundaries(domain))))
    w = vertex(domain, mod1(i + 1, length(boundaries(domain))))

    if distribution == :linear
        f = j -> domain.parent(j // (n + 1))
    elseif distribution == :chebyshev
        f = j -> (cospi(domain.parent((2j - 1) // 2n)) + 1) / 2
    elseif distribution == :exponential
        m = (n + 1) // 2
        f = j -> begin
            d = exp(-domain.parent(m - j))
            if j < m
                return d / 2
            elseif j > m
                return 1 - inv(d) / 2
            else
                return domain.parent(0.5)
            end
        end
    elseif distribution == :root_exponential
        m = (n + 1) // 2
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
        throw(ArgumentError("no distribution named $distribution"))
    end

    points = Vector{SVector{2,arb}}(undef, n)
    for j = 1:n
        t = f(j)
        points[j] = v .+ t .* (w - v)
    end

    return points, fill(i, n)
end

# TODO: Fix this
function interior_points(domain::Polygon, n::Integer; rng = MersenneTwister(42))
    points = Vector{SVector{2,arb}}(undef, n)
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
