function Triangle{T}(domain::Triangle{T}, parent::ArbField) where {T <: Union{arb, fmpq}}
    Triangle(domain.α, domain.β, parent)
end

function Base.getproperty(domain::Triangle{T}, name::Symbol) where {T}
    if name == :γ
        if T == fmpq
            return 1 - domain.α - domain.β
        else
            return domain.parent(π) - domain.α - domain.β
        end
    else
        return getfield(domain, name)
    end
end

function Base.show(io::IO, domain::Triangle{fmpq})
    angles_string = ["$(numerator(a))π/$(denominator(a))" for a in (domain.α, domain.β, domain.γ)]
    print(io, "Triangle with $(domain.parent.prec) bits of precision and angles ($(angles_string[1]), $(angles_string[2]), $(angles_string[3]))")
end

function Base.show(io::IO, domain::Triangle{arb})
    angles_string = ["$a" for a in (domain.α, domain.β, domain.γ)]
    print(io, "Triangle with $(domain.parent.prec) bits of precision and angles ($(angles_string[1]), $(angles_string[2]), $(angles_string[3]))")
end

boundaries(::Triangle) = 1:3

"""
    angle(domain::Triangle, i::Integer)

Return the angle for vertex `i` of the triangle.
"""
function angle(domain::Triangle{fmpq}, i::Integer)
    if i == 1
        return domain.parent(π)*domain.α
    elseif i == 2
        return domain.parent(π)*domain.β
    elseif i == 3
        return domain.parent(π)*domain.γ
    else
        throw(ArgumentError("attempt to get vertex $i from a $(typeof(domain))"))
    end
end

function angle(domain::Triangle{arb}, i::Integer)
    if i == 1
        return domain.α
    elseif i == 2
        return domain.β
    elseif i == 3
        return domain.γ
    else
        throw(ArgumentError("attempt to get vertex $i from a $(typeof(domain))"))
    end
end

"""
    angledivπ(domain::Triangle, i::Integer)

Return the angle for vertex `i` of the triangle divided by `π`.
"""
function angledivπ(domain::Triangle{fmpq}, i::Integer)
    if i == 1
        return domain.α
    elseif i == 2
        return domain.β
    elseif i == 3
        return domain.γ
    else
        throw(ArgumentError("attempt to get vertex $i from a $(typeof(domain))"))
    end
end

function angledivπ(domain::Triangle{arb}, i::Integer)
    if i == 1
        return domain.α/domain.parent(π)
    elseif i == 2
        return domain.β/domain.parent(π)
    elseif i == 3
        return domain.γ/domain.parent(π)
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
            x = sin(angle(domain, 2))/sin(angle(domain, 3))
            s, c = sincos(angle(domain, 1))
        else
            x = sinpi(domain.β, domain.parent)/sinpi(domain.γ, domain.parent)
            s, c = sincospi(domain.α, domain.parent)
        end
        return SVector(c*x, s*x)
    else
        throw(ArgumentError("attempt to get vertex $i from a $(typeof(domain))"))
    end
end

"""
    vertices(domain::Triangle)

Return all three vertices of the triangle.
"""
vertices(domain::Triangle) = tuple((vertex(domain, i) for i in boundaries(domain))...)

area(domain::Triangle) = vertex(domain, 3)[2]/2

"""
    center(domain::Triangle)

Return the centroid of the triangle.
"""
center(domain::Triangle) = sum(vertices(domain))/3

function boundary_parameterization(t, domain::Triangle, i::Integer)
    i ∈ boundaries(domain) || throw(ArgumentError("attempt to get vertex $i from a $(typeof(domain))"))
    v = vertex(domain, mod1(i + 1, 3))
    w = vertex(domain, mod1(i + 2, 3))
    return v .+ t.*(w - v)
end

function boundary_points(domain::Triangle, i::Integer, n::Integer)
    i ∈ boundaries(domain) || throw(ArgumentError("attempt to get vertex $i from a $(typeof(domain))"))

    v = vertex(domain, mod1(i + 1, 3))
    w = vertex(domain, mod1(i + 2, 3))

    points = Vector{SVector{2, arb}}(undef, n)
    for j in 1:n
        t = domain.parent(j//(n + 1))
        points[j] = v .+ t.*(w - v)
    end

    return points, fill(i, n)
end

function interior_points(domain::Triangle, n::Integer; rng = MersenneTwister(42))
    v = vertex(domain, 2)
    w = vertex(domain, 3)

    points = Vector{SVector{2, arb}}(undef, n)
    for i in 1:n
        u1, u2 = rand(rng), rand(rng)
        if u1 + u2 > 1
            u1 = 1 - u1
            u2 = 1 - u2
        end
        points[i] = u1.*v + u2.*w
    end

    return points
end
