function Triangle{S,T}(domain::Triangle{S,T}; parent::ArbField) where {S,T}
    Triangle(domain.angles[1], domain.angles[2]; parent)
end

function Base.show(io::IO, domain::Triangle{S,T}) where {S,T}
    if has_rational_angles(domain)
        angles_string = ["$(numerator(a))π/$(denominator(a))" for a in domain.angles]
    else
        angles_string = ["$a" for a in domain.angles]
    end
    print(
        io,
        "Triangle{$S,$T} with  angles ($(angles_string[1]), $(angles_string[2]), $(angles_string[3]))",
    )
    if !isnothing(domain.parent)
        print(io, " and $(precision(domain.parent)) bits of precision")
    end
end

vertexindices(::Triangle) = 1:3
boundaries(::Triangle) = 1:3

angle_raw(domain::Triangle, i::Integer) = domain.angles[i]

function vertex(domain::Triangle{S,T}, i::Integer) where {S,T}
    if i == 1
        return SVector(zero(S), zero(S))
    elseif i == 2
        return SVector(one(S), zero(S))
    elseif i == 3
        if has_rational_angles(domain)
            x =
                sinpi(convert(S, angle_raw(domain, 2))) /
                sinpi(convert(S, angle_raw(domain, 3)))
            s, c = sincospi(convert(S, angle_raw(domain, 1)))
        else
            x = sin(angle(domain, 2)) / sin(angle(domain, 3))
            s, c = sincos(angle(domain, 1))
        end
        return SVector(c * x, s * x)
    else
        throw(ArgumentError("attempt to get vertex $i from a $(typeof(domain))"))
    end
end

function vertex(domain::Triangle{arb,T}, i::Integer) where {T}
    if i == 1
        return SVector(domain.parent(0), domain.parent(0))
    elseif i == 2
        return SVector(domain.parent(1), domain.parent(0))
    elseif i == 3
        if has_rational_angles(domain)
            x =
                sinpi(angle_raw(domain, 2), domain.parent) /
                sinpi(angle_raw(domain, 3), domain.parent)
            s, c = sincospi(angle_raw(domain, 1), domain.parent)
        else
            x = sin(angle(domain, 2)) / sin(angle(domain, 3))
            s, c = sincos(angle(domain, 1))
        end
        return SVector(c * x, s * x)
    else
        throw(ArgumentError("attempt to get vertex $i from a $(typeof(domain))"))
    end
end

function orientation_raw(domain::Triangle{S,T}, i::Integer; reversed = false) where {S,T}
    π = has_rational_angles(domain) ? one(T) : pi

    if i == 1
        res = zero(T)
    elseif i == 2
        res = π - angle_raw(domain, 2)
    elseif i == 3
        res = π + angle_raw(domain, 1)
    else
        throw(ArgumentError("attempt to get vertex $i from a $(typeof(domain))"))
    end

    if reversed
        return 2convert(T, π) - angle_raw(domain, i) - res
    else
        return res
    end
end

function orientation_raw(domain::Triangle{arb,T}, i::Integer; reversed = false) where {T}
    if T == fmpq
        π = fmpq(1)
    else
        π = domain.parent(pi)
    end

    if i == 1
        res = T == fmpq ? fmpq(0) : domain.parent(0)
    elseif i == 2
        res = π - angle_raw(domain, 2)
    elseif i == 3
        res = π + angle_raw(domain, 1)
    else
        throw(ArgumentError("attempt to get vertex $i from a $(typeof(domain))"))
    end

    if reversed
        return 2π - angle_raw(domain, i) - res
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
