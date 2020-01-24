function Base.show(io::IO, domain::SphericalTriangle{fmpq})
    angles_string = ["$(numerator(a))π/$(denominator(a))" for a in domain.angles]
    print(io, "Spherical triangle with angles ($(angles_string[1]), $(angles_string[2]), $(angles_string[3]))")
end

function Base.show(io::IO, domain::SphericalTriangle{arb})
    print(io, "Spherical triangle with angles $(domain.angles)")
end

"""
    angles(domain::SphericalTriangle)
> Return the angles of the spherical triangle.
"""
function angles(domain::SphericalTriangle{fmpq})
    domain.parent(π).*domain.angles
end

function angles(domain::SphericalTriangle{arb})
    domain.angles
end

"""
    angle(domain::SphericalTriangle, i::Integer)
> Return the angle for vertex i of the spherical triangle.
"""
function angle(domain::SphericalTriangle{fmpq}, i::Integer)
    domain.parent(π)*domain.angles[i]
end

function angle(domain::SphericalTriangle{arb}, i::Integer)
    domain.angles[i]
end

"""
    vertex(domain::SphericalTriangle, i::Integer)
> Return cartesian coordinates for vertex i of the spherical triangle.
"""
function vertex(domain::SphericalTriangle, i::Integer)
    RR = domain.parent
    if i == 1
        return RR.(SVector(0, 0, 1))
    elseif i == 2
        α, β, γ = angles(domain)
        S = RR(α + β + γ)/2
        θ = 2asin(sqrt(-cos(S)*cos(S - γ)/(sin(RR(α))*sin(RR(β)))))
        return SVector(sin(θ), RR(0), cos(θ))
    elseif i == 3
        α, β, γ = angles(domain)
        S = RR(α + β + γ)/2
        θ = 2asin(sqrt(-cos(S)*cos(S - β)/(sin(RR(α))*sin(RR(γ)))))
        return SVector(sin(θ)*cos(RR(α)), sin(θ)*sin(RR(α)), cos(θ))
    else
        throw(ErrorException("attempt to get vertex number $i from a spherical triangle"))
    end
end

"""
    area(domain::SphericalTriangle)
> Compute the area of the spherical triangle.
"""
function area(domain::SphericalTriangle)
    sum(angles(domain)) - domain.parent(π)
end

"""
    center(domain::SphericalTriangle)
> Return (θ, ϕ) for the center of the spherical triangle, given by the
  normalised sum of the vertices.
"""
function center(domain::SphericalTriangle)
    spherical(sum(vertex(domain, i) for i in 1:3))
end

"""
    greatcircleplane(domain::SphericalTriangle)
> Return coefficients a, b, c giving the plane ax + by + cz = 0 which
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
Compute the θ value corresponding to the give ϕ value on the great
circle defined by the plane ax + by + cz = 0.
"""
function greatcircle(ϕ::arb, a::arb, b::arb, c::arb)
    θ = atan(-c, a*cos(ϕ) + b*sin(ϕ))
    if θ < 0
        θ += ϕ.parent(π)
    end
    return θ
end

"""
    boundary_points(domain::SphericalTriangle,
                    i::Integer,
                    n::Integer)
> Return (θ, ϕ) for n points on the boundary opposite of vertex number
  i.
"""
function boundary_points(domain::SphericalTriangle,
                         i::Integer,
                         n::Integer)
    i >= 1 && i <= 3 || throw(ErrorException("attempt to use vertex number $i from a spherical triangle"))
    RR = domain.parent

    points = Vector{NamedTuple{(:θ, :ϕ),Tuple{arb,arb}}}(undef, n)

    v = vertex(domain, mod1(i + 1, 3))
    w = vertex(domain, mod1(i + 2, 3))

    for j in 1:n
        t = RR(j//(n + 1))

        xyz = v .+ t.*(w - v)
        points[j] = spherical(xyz)
    end

    points
end

"""
    boundary_points(domain::SphericalTriangle,
                    eigenfunction::AbstractSphericalEigenfunction,
                    n::Integer)
> Return (θ, ϕ) for n/3 points on each boundary of the spherical
  triangle. They are returned in order, so first n/3 points from the
  first boundary, then n/3 from the second and finally n/3 from the
  third. In case n is not divisible by 3 then it takes one extra point
  from the first and possibly the second boundary.
"""
function boundary_points(domain::SphericalTriangle,
                         eigenfunction::AbstractSphericalEigenfunction,
                         n::Integer)
    [boundary_points(domain, 1, div(n, 3) + ifelse(n % 3 >= 1, 1, 0));
     boundary_points(domain, 2, div(n, 3) + ifelse(n % 3 >= 2, 1, 0));
     boundary_points(domain, 3, div(n, 3))]
end

"""
    boundary_points(domain::SphericalTriangle,
                    eigenfunction::SphericalVertexEigenfunction,
                    n::Integer)
> Return (θ, ϕ) for n points on the boundary opposite of the vertex
  used for the eigenfunction.
"""
function boundary_points(domain::SphericalTriangle,
                         eigenfunction::SphericalVertexEigenfunction,
                         n::Integer)
    boundary_points(domain, eigenfunction.vertex, n)
end

"""
    interior_points(domain::SphericalTriangle, n::Integer; rng = MersenneTwister(42))
> Return (θ, ϕ) for n points taken uniformly at random from the
  interior of the triangle.
"""
function interior_points(domain::SphericalTriangle,
                         n::Integer;
                         rng = MersenneTwister(42))
    RR = domain.parent

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

    points = Vector{NamedTuple{(:θ, :ϕ),Tuple{arb,arb}}}(undef, n)

    # Maximum value of ϕ
    A = angle(domain, 1)

    # Compute coefficients for plane determining the great circle
    # giving the lower boundary of the spherical triangle.
    (d, e, f) = greatcircleplane(domain)

    for i in 1:n
        # Randomly pick θ and ϕ
        θ = acos(1 - 2RR(rand(rng)))
        ϕ = rand(rng)*A

        # To determine if the point (θ, ϕ) is inside the triangle we
        # check if θ is above the θ value of great circle going
        # between the two vertices of the triangle not on the north
        # pole at the value ϕ.
        while θ > greatcircle(ϕ, d, e, f)
            θ = acos(1 - 2RR(rand(rng)))
            ϕ = rand(rng)*A
        end

        points[i] = (θ = θ, ϕ = ϕ)
    end

    points
end
