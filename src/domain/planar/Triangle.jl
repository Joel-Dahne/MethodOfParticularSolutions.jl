function Triangle{T}(domain::Triangle{T}, parent::ArbField) where {T<:Union{arb,fmpq}}
    Triangle(domain.angles[1], domain.angles[2], parent)
end

function Base.show(io::IO, domain::Triangle{fmpq})
    angles_string = ["$(numerator(a))π/$(denominator(a))" for a in domain.angles]
    print(
        io,
        "Triangle with $(domain.parent.prec) bits of precision and angles ($(angles_string[1]), $(angles_string[2]), $(angles_string[3]))",
    )
end

function Base.show(io::IO, domain::Triangle{arb})
    angles_string = ["$a" for a in domain.angles]
    print(
        io,
        "Triangle with $(domain.parent.prec) bits of precision and angles ($(angles_string[1]), $(angles_string[2]), $(angles_string[3]))",
    )
end

vertexindices(::Triangle) = 1:3
boundaries(::Triangle) = 1:3

angle_raw(domain::Triangle, i::Integer) = domain.angles[i]

"""
    angle(domain::Triangle, i::Integer)

Return the angle for vertex `i` of the triangle.
"""
function angle(domain::Triangle{fmpq}, i::Integer)
    if i ∈ boundaries(domain)
        return domain.parent(π) * domain.angles[i]
    else
        throw(ArgumentError("attempt to get vertex $i from a $(typeof(domain))"))
    end
end

function angle(domain::Triangle{arb}, i::Integer)
    if i ∈ boundaries(domain)
        return domain.angles[i]
    else
        throw(ArgumentError("attempt to get vertex $i from a $(typeof(domain))"))
    end
end

"""
    angledivπ(domain::Triangle, i::Integer)

Return the angle for vertex `i` of the triangle divided by `π`.
"""
function angledivπ(domain::Triangle{fmpq}, i::Integer)
    if i ∈ boundaries(domain)
        return domain.angles[i]
    else
        throw(ArgumentError("attempt to get vertex $i from a $(typeof(domain))"))
    end
end

function angledivπ(domain::Triangle{arb}, i::Integer)
    if i ∈ boundaries(domain)
        return domain.angles[i] / domain.parent(π)
    else
        throw(ArgumentError("attempt to get vertex $i from a $(typeof(domain))"))
    end
end

"""
    angles(domain::Triangle)

Return the angles of the triangle.
"""
angles(domain::Triangle) = tuple((angle(domain, i) for i in boundaries(domain))...)

"""
    anglesdivπ(domain::Triangle)

Return the angles of the triangle divided by `π`.
"""
anglesdivπ(domain::Triangle) = tuple((angledivπ(domain, i) for i in boundaries(domain))...)

"""
    vertex(domain::Triangle, i::Integer)

Return the Cartesian coordinates of vertex `i` of the triangle.
"""
function vertex(domain::Triangle{T}, i::Integer) where {T}
    if i == 1
        return SVector(domain.parent(0), domain.parent(0))
    elseif i == 2
        return SVector(domain.parent(1), domain.parent(0))
    elseif i == 3
        if T == arb
            x = sin(angle(domain, 2)) / sin(angle(domain, 3))
            s, c = sincos(angle(domain, 1))
        else
            x =
                sinpi(angledivπ(domain, 2), domain.parent) /
                sinpi(angledivπ(domain, 3), domain.parent)
            s, c = sincospi(angledivπ(domain, 1), domain.parent)
        end
        return SVector(c * x, s * x)
    else
        throw(ArgumentError("attempt to get vertex $i from a $(typeof(domain))"))
    end
end

"""
    vertices(domain::Triangle)

Return all three vertices of the triangle.
"""
vertices(domain::Triangle) = tuple((vertex(domain, i) for i in boundaries(domain))...)

function orientation(domain::Triangle{T}, i::Integer; reversed = false) where {T}
    if T == fmpq
        π = fmpq(1)
    else
        π = domain.parent(pi)
    end

    if i == 1
        res = T == fmpq ? fmpq(0) : domain.parent(0)
    elseif i == 2
        res = π - domain.angles[2]
    elseif i == 3
        res = π + domain.angles[1]
    else
        throw(ArgumentError("attempt to get vertex $i from a $(typeof(domain))"))
    end

    if reversed
        return 2π - domain.angles[i] - res
    else
        return res
    end
end

area(domain::Triangle) = vertex(domain, 3)[2] / 2

"""
    center(domain::Triangle)

Return the centroid of the triangle.
"""
center(domain::Triangle) = sum(vertices(domain)) / 3

function Base.in(xy, domain::Triangle)
    x, y = xy
    return (
        y >= 0 &&
        atan(y, x) <= angle(domain, 1) &&
        domain.parent(π) - atan(y, x - 1) <= angle(domain, 2) # To the left of the right edge
    )
end

function boundary_parameterization(t, domain::Triangle, i::Integer)
    i ∈ boundaries(domain) ||
        throw(ArgumentError("attempt to get vertex $i from a $(typeof(domain))"))
    v = vertex(domain, mod1(i + 1, 3))
    w = vertex(domain, mod1(i + 2, 3))
    return v .+ t .* (w - v)
end

function boundary_points(
    domain::Triangle,
    i::Integer,
    n::Integer;
    distribution = :chebyshev,
)
    i ∈ boundaries(domain) ||
        throw(ArgumentError("attempt to get vertex $i from a $(typeof(domain))"))

    v = vertex(domain, mod1(i + 1, 3))
    w = vertex(domain, mod1(i + 2, 3))

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
        f = j -> begin
            d = exp(-domain.parent(m - j) / sqrt(domain.parent(n)))
            if j < m
                return d / 2
            elseif j > m
                return 1 - inv(d) / 2
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

function interior_points(domain::Triangle, n::Integer; rng = MersenneTwister(42))
    v = vertex(domain, 2)
    w = vertex(domain, 3)

    points = Vector{SVector{2,arb}}(undef, n)
    for i = 1:n
        u1, u2 = rand(rng), rand(rng)
        if u1 + u2 > 1
            u1 = 1 - u1
            u2 = 1 - u2
        end
        points[i] = u1 .* v + u2 .* w
    end

    return points
end
