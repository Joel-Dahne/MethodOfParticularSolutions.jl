function SphericalTriangle{T}(
    domain::SphericalTriangle;
    parent::ArbField,
) where {T<:Union{fmpq,arb}}
    SphericalTriangle(domain.angles; parent)
end

function Base.show(io::IO, domain::SphericalTriangle{fmpq})
    angles_string = ["$(numerator(a))π/$(denominator(a))" for a in domain.angles]
    print(
        io,
        "Spherical triangle with angles ($(angles_string[1]), $(angles_string[2]), $(angles_string[3]))",
    )
end

function Base.show(io::IO, domain::SphericalTriangle{arb})
    print(io, "Spherical triangle with angles $(domain.angles)")
end

"""
    angles(domain::SphericalTriangle)
> Return the angles of the spherical triangle.
"""
function angles(domain::SphericalTriangle{fmpq})
    domain.parent(π) .* domain.angles
end

function angles(domain::SphericalTriangle{arb})
    domain.angles
end

"""
    angle(domain::SphericalTriangle, i::Integer)
> Return the angle for vertex i of the spherical triangle.
"""
function angle(domain::SphericalTriangle{fmpq}, i::Integer)
    domain.parent(π) * domain.angles[i]
end

function angle(domain::SphericalTriangle{arb}, i::Integer)
    domain.angles[i]
end

"""
    vertex(domain::SphericalTriangle, i::Integer)
> Compute vertex i of the spherical triangle.
"""
function vertex(domain::SphericalTriangle, i::Integer)
    RR = domain.parent
    if i == 1
        return RR.(SVector(0, 0, 1))
    elseif i == 2
        α, β, γ = angles(domain)
        S = (α + β + γ) / 2
        θ = 2asin(sqrt(-cos(S) * cos(S - γ) / (sin(α) * sin(β))))
        return SVector(sin(θ), RR(0), cos(θ))
    elseif i == 3
        α, β, γ = angles(domain)
        S = (α + β + γ) / 2
        θ = 2asin(sqrt(-cos(S) * cos(S - β) / (sin(α) * sin(γ))))
        return SVector(sin(θ) * cos(α), sin(θ) * sin(α), cos(θ))
    else
        throw(ErrorException("attempt to get vertex $i from a spherical triangle"))
    end
end

function area(domain::SphericalTriangle)
    sum(angles(domain)) - domain.parent(π)
end

"""
    center(domain::SphericalTriangle)
> Compute the center of the spherical triangle, given by the
  normalised sum of the vertices.
"""
function center(domain::SphericalTriangle)
    xyz = sum(vertex(domain, i) for i = 1:3)
    r = sqrt(sum(xyz .^ 2))
    xyz ./ r
end

"""
    greatcircleplane(domain::SphericalTriangle)
> Compute coefficients a, b, c giving the plane ax + by + cz = 0 which
  is parallel to the great circle that defines the bottom boundary of
  the spherical triangle.
"""
function greatcircleplane(domain::SphericalTriangle)
    # Compute the normal vector of the plane
    w = cross(vertex(domain, 2), vertex(domain, 3))
    # w = [a, b, c]
    return (w[1], w[2], w[3])
end

"""
    greatcircle(ϕ::arb, a::arb, b::arb, c::arb)
> Compute the θ value corresponding to the give ϕ value on the great
circle defined by the plane ax + by + cz = 0.
"""
function greatcircle(ϕ::arb, a::arb, b::arb, c::arb)
    θ = atan(-c, a * cos(ϕ) + b * sin(ϕ))
    if θ < 0
        θ += ϕ.parent(π)
    end
    return θ
end

"""
    theta_bound(domain::SphericalTriangle, type = :lower)
> Compute a bound of θ on the boundary opposite of the north pole. By
  default a lower bound is computed, an upper bound is computed with
  `type = :upper`.
"""
function theta_bound(domain::SphericalTriangle, type = :lower)
    if type == :lower
        cmp = min
    elseif type == :upper
        cmp = max
    else
        throw(ArgumentError("type must be :lower or :upper, got $type"))
    end

    # Compute the critical point with respect to θ on the boundary
    i = 1
    u = vertex(domain, 2)
    v = vertex(domain, 3)
    w = u .- v

    critical_point = (
        (v[3] * sum(v .* w) - w[3] * sum(v .* v)) /
        (w[3] * sum(v .* w) - v[3] * sum(w .* w))
    )

    if !isfinite(critical_point)
        @warn "Problem computing critical point in theta_bound, got $critical_point"
    end

    a, b, c = greatcircleplane(domain)
    # Compute θ for the two endpoints
    θ = cmp(greatcircle(spherical(u).ϕ, a, b, c), greatcircle(spherical(v).ϕ, a, b, c))

    # If the critical point lies in the interval (0, 1) compute the
    # angle corresponding to it
    if !isnonpositive(critical_point) && !isnonnegative(critical_point - 1)
        θ = cmp(θ, greatcircle(spherical(v .+ critical_point .* w).ϕ, a, b, c))
    end

    return θ
end

function boundary_parameterization(t, domain::SphericalTriangle, i::Integer)
    i >= 1 && i <= 3 ||
        throw(ErrorException("attempt to use vertex number $i from a spherical triangle"))
    v = vertex(domain, mod1(i + 1, 3))
    w = vertex(domain, mod1(i + 2, 3))

    normalize(v .+ t .* (w - v))
end

"""
    boundary_points(domain::SphericalTriangle,
                    i::Integer,
                    n::Integer;
                    distribution = :linear,
                    half::Bool = true)
> Return n points on the boundary opposite of vertex number i.

  If `half` is true then only take points from the first half of the
  boundary.
"""
function boundary_points(
    domain::SphericalTriangle,
    i::Integer,
    n::Integer;
    distribution = :linear,
    half::Bool = false,
)
    i >= 1 && i <= 3 ||
        throw(ErrorException("attempt to use vertex number $i from a spherical triangle"))
    points = Vector{SVector{3,arb}}(undef, n)

    v = vertex(domain, mod1(i + 1, 3))
    w = vertex(domain, mod1(i + 2, 3))

    for j = 1:n
        t = ifelse(half, domain.parent(j // (2n + 1)), domain.parent(j // (n + 1)))

        points[j] = normalize(v .+ t .* (w - v))
    end

    return points, fill(i, n)
end

"""
    boundary_points(domain::SphericalTriangle,
                    eigenfunction::SphericalVertexEigenfunction,
                    n::Integer)
> Return n points on the boundary opposite of the vertex used for the
  eigenfunction.
"""
function boundary_points(
    domain::SphericalTriangle,
    eigenfunction::SphericalVertexEigenfunction,
    n::Integer,
)
    boundary_points(domain, eigenfunction.vertex, n, half = eigenfunction.stride == 2)
end

"""
    boundary_points(domain::SphericalTriangle,
                    eigenfunction::KrewerasEigenfunction,
                    n::Integer)
> Return n points from half of boundary one.
"""
function boundary_points(
    domain::SphericalTriangle,
    eigenfunction::KrewerasEigenfunction,
    n::Integer,
)
    boundary_points(domain, 1, n, half = true)
end

function interior_points(domain::SphericalTriangle, n::Integer; rng = MersenneTwister(42))
    # The points are generated by uniformly choosing a value for ϕ and
    # for θ we uniformly choose an u ∈ [-1, 1] and let θ = acos(1 -
    # u), this gives a uniform distribution on the sphere. To ensure
    # that the generated point is indeed inside the triangle we
    # restrict the possible values for ϕ, given by the angle at the
    # north pole. For θ we generate it and then check that it's inside
    # the triangle by checking that it's above the great circle
    # determining the lower boundary of the spherical triangle, if
    # it's outside we generate ϕ and θ again until we get something
    # which is inside.
    RR = domain.parent

    points = Vector{SVector{3,arb}}(undef, n)

    # Maximum value of ϕ
    A = angle(domain, 1)

    # Compute coefficients for plane determining the great circle
    # giving the lower boundary of the spherical triangle.
    (d, e, f) = greatcircleplane(domain)

    for i = 1:n
        # Randomly pick θ and ϕ
        θ = acos(1 - 2RR(rand(rng)))
        ϕ = rand(rng) * A

        # To determine if the point (θ, ϕ) is inside the triangle we
        # check if θ is above the θ value of great circle going
        # between the two vertices of the triangle not on the north
        # pole at the value ϕ.
        while θ > greatcircle(ϕ, d, e, f)
            θ = acos(1 - 2RR(rand(rng)))
            ϕ = rand(rng) * A
        end

        points[i] = cartesian(θ, ϕ)
    end

    points
end

"""
    anglesfromvertices(a, b, c)
> Compute the angles for the spherical triangle having the vertices
  `a`, `b` and `c`.
"""
function anglesfromvertices(a, b, c)
    v = cross(a, b)
    w = cross(a, c)
    α = acos(dot(v, w) / LinearAlgebra.norm(v) / LinearAlgebra.norm(w))

    v = cross(b, a)
    w = cross(b, c)
    β = acos(dot(v, w) / LinearAlgebra.norm(v) / LinearAlgebra.norm(w))

    v = cross(c, a)
    w = cross(c, b)
    γ = acos(dot(v, w) / LinearAlgebra.norm(v) / LinearAlgebra.norm(w))

    (α, β, γ)
end

"""
    subtriangle(domain::SphericalTriangle, ratio = 0.5)
> Return vertices for a spherical triangle in the interior of the
  domain. The size of the triangle is determined by `ratio`.
"""
function subtriangle(domain::SphericalTriangle; ratio = 0.5)
    v = center(domain)

    a = normalize((1 - ratio) * v + ratio * vertex(domain, 1))
    b = normalize((1 - ratio) * v + ratio * vertex(domain, 2))
    c = normalize((1 - ratio) * v + ratio * vertex(domain, 3))

    return (a, b, c)
end

"""
    partitiontriangle(a, b, c; iterations::Int = 1)
> Given the vertices of a spherical triangle return vertices for the
  four spherical triangle obtained by partitioning the original
  triangle into smaller ones.

  If `iterations` is not equal to 1 then recursively do this the give
  number of times, in case `iterations = 0` then return the original
  triangle.
"""
function partitiontriangle(a, b, c; iterations::Int = 1)
    if iterations < 1
        return [(a, b, c)]
    end
    x = normalize(a + b)
    y = normalize(b + c)
    z = normalize(c + a)

    triangles = [(a, x, z), (x, b, y), (z, y, c), (x, y, z)]
    if iterations == 1
        return triangles
    else
        return collect(
            Iterators.flatten([
                partitiontriangle(i, j, k, iterations = iterations - 1) for
                (i, j, k) in triangles
            ]),
        )
    end
end
