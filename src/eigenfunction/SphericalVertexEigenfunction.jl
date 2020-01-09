function SphericalVertexEigenfunction(domain::SphericalTriangle,
                                      vertex::Int)
    vertex >= 1 && vertex <= 3 || throw(ErrorException("vertex must be between 1 and 3 not $vertex"))
    SphericalVertexEigenfunction(domain, vertex, arb[])
end

"""
    coordinate_transformation(u::SphericalVertexEigenfunction)
> Return a coordinate transformation T which switches from spherical
  coordinates (θ, ϕ) to (θ', ϕ'), in the new coordinate system the
  vertex from which the eigenfunction originates is put at the north
  pole with one of the adjacent boundaries being parallel to the
  y-axis.

  This is computed by switching to cartesian coordinates, performing a
  number of rotations and then switching back to spherical
  coordinates.
"""
function coordinate_transformation(u::SphericalVertexEigenfunction)
    # TODO: The performance could most likely be improved
    if u.vertex == 1
        T = (θ, ϕ) -> begin
            θ, ϕ
        end
    elseif u.vertex == 2
        # Rotate by α along the y-axis so that the second vertex ends
        # up at the north pole, then rotate by β around the z-axis so
        # that the boundary ends up parallel to the y-axis.
        α = -acos(vertex(u.domain, 2)[3])
        β = angles(u.domain)[2] - u.domain.parent(π)
        L = LinearMap(RotZY(β, α))

        T = (θ, ϕ) -> begin
            xyz = cartesian(θ, ϕ)
            xyz = L(xyz)
            (θ, ϕ) = spherical(xyz)
            return θ, ϕ
        end
    elseif u.vertex == 3
        # Rotate by α along the z-axis so that the third vertex
        # ends up parallel to the y-axis, then rotate by β around
        # the y-axis so that it ends up at the north pole. Finally
        # rotate by γ so along the z-axis so that the boundary
        # ends up parallel to the y-axis.
        α = -angles(u.domain)[1]
        β = -acos(vertex(u.domain, 3)[3])
        γ = u.domain.parent(π)
        L = LinearMap(RotZYZ(γ, β, α))

        T = (θ, ϕ) -> begin
            xyz = cartesian(θ, ϕ)
            xyz = L(xyz)
            (θ, ϕ) = spherical(xyz)
            return θ, ϕ
        end
    else
        throw(ErrorException("vertex must be between 1 and 3 not $vertex"))
    end

    return T
end

"""
    mu(eigenfunction::SphericalVertexEigenfunction,
       k::Integer = 1)
> Return k*μ0 as an arb ball, which is the parameter used for the
  Legendre function.
"""
function mu(u::SphericalVertexEigenfunction{fmpq},
            k::Integer = 1)
    u.domain.parent(-k*inv(u.domain.angles[u.vertex]))
end

function mu(u::SphericalVertexEigenfunction{arb},
            k::Integer = 1)
    u.domain.parent(-k*inv(u.domain.angles[u.vertex]))*u.domain.parent(π)
end

"""
    u(θ::arb, ϕ::arb, λ::arb, k::Integer; notransform::Bool = false)
> Evaluate the k-th basis function for the eigenfunction with the
  given λ on the point given by (θ, ϕ). If notransform is true then
  do not perform a coordinate transform on θ and ϕ first, this assumes
  that they already given in the coordinate system used by u. See
  coordinate_transform for details about the coordinate transform
  used.
"""
function (u::SphericalVertexEigenfunction)(θ::arb,
                                           ϕ::arb,
                                           λ::arb,
                                           k::Integer;
                                           notransform::Bool = false)
    ν = θ.parent(-0.5) + sqrt(θ.parent(0.25) + λ)
    μ = θ.parent(mu(u, k))
    if !notransform
        θ, ϕ = coordinate_transformation(u, θ, ϕ)
    end
    legendre_p_safe(ν, μ, cos(θ))*sin(μ*ϕ)
end

"""
    u(θ::arb, ϕ::arb, λ::arb; notransform::Bool = false)
> Evaluate the eigenfunction with the given λ on the point given by
  (θ, ϕ). If notransform is true then do not perform a coordinate
  transform on θ and ϕ first, this assumes that they already given in
  the coordinate system used by u. See coordinate_transform for
  details about the coordinate transform used.
"""
function (u::SphericalVertexEigenfunction)(θ::arb,
                                           ϕ::arb,
                                           λ::arb;
                                           notransform::Bool = false)
    res = θ.parent(0)

    if !notransform
        θ, ϕ = coordinate_transformation(u, θ, ϕ)
    end

    for k in 1:length(u.coefficients)
        res += u.coefficients[k]*u(θ, ϕ, λ, k, notransform = true)
    end

    res
end
