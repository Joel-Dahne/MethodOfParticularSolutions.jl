function Polygon{T}(domain::Polygon{T}, parent::ArbField = domain.parent) where {T <: Union{arb, fmpq}}
    return Polygon(copy(domain.angles), copy(domain.vertices), parent)
end


function Polygon(angles::AbstractVector{T}, vertices::AbstractVector, parent::ArbField) where {T}
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
    print(io, "Polygon with $(length(boundaries(domain))) sides and $(domain.parent.prec) bits of precision")
end

boundaries(domain::Polygon) = 1:length(domain.angles)

angle(domain::Polygon{fmpq}, i::Integer) = domain.parent(π)*domain.angles[i]
angle(domain::Polygon{arb}, i::Integer) = domain.angles[i]

angledivπ(domain::Polygon{fmpq}, i::Integer) = domain.angles[i]
angledivπ(domain::Polygon{arb}, i::Integer) = domain.angles[i]/domain.parent(π)

angles(domain::Polygon) = [angle(domain, i) for i in boundaries(domain)]
anglesdivπ(domain::Polygon) = [angledivπ(domain, i) for i in boundaries(domain)]

vertex(domain::Polygon, i::Integer) = domain.vertices[i]
vertices(domain::Polygon) = domain.vertices

# TODO: Fix this
area(domain::Polygon) = 1
center(domain::Polygon) = sum(vertices(domain))/length(boundaries(domain))

# TODO: Add this
function Base.in(xy, domain::Polygon)
    vs = vertices(domain)

    isleft(p0, p1, p2) = (p1[1] - p0[1])*(p2[2] - p0[2]) - (p2[1] - p0[1])*(p1[2] - p0[2])

    wn = 0
    for i in eachindex(vs)
        vi = vs[i]
        vip1 = vs[mod1(i + 1, length(vs))]
        if vi[2] <= xy[2]
            if vip1[2] > xy[2]
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
    i ∈ boundaries(domain) || throw(ArgumentError("attempt to get vertex $i from a $(typeof(domain))"))
    v = vertex(domain, mod1(i, length(boundaries(domain))))
    w = vertex(domain, mod1(i + 1, length(boundaries(domain))))

    return v .+ t.*(w - v)
end

function boundary_points(domain::Polygon, i::Integer, n::Integer)
    i ∈ boundaries(domain) || throw(ArgumentError("attempt to get vertex $i from a $(typeof(domain))"))

    v = vertex(domain, mod1(i, length(boundaries(domain))))
    w = vertex(domain, mod1(i + 1, length(boundaries(domain))))

    points = Vector{SVector{2, arb}}(undef, n)
    for j in 1:n
        t = domain.parent(j//(n + 1))
        points[j] = v .+ t.*(w - v)
    end

    return points, fill(i, n)
end

# TODO: Fix this
function interior_points(domain::Polygon, n::Integer; rng = MersenneTwister(42))
    points = Vector{SVector{2, arb}}(undef, n)
    xmin, xmax = extrema(getindex.(vertices(domain), 1))
    ymin, ymax = extrema(getindex.(vertices(domain), 2))
    for i in 1:n
        pt = SVector(xmin + rand(rng)*(xmax - xmin), ymin + rand(rng)*(ymax - ymin))
        while !(pt ∈ domain)
            pt = SVector(xmin + rand(rng)*(xmax - xmin), ymin + rand(rng)*(ymax - ymin))
        end

        points[i] = pt
    end

    return points
end