function active_boundaries(u::AbstractSphericalEigenfunction)
    (true, true, true)
end

"""
    coordinate_transformation(u::AbstractSphericalEigenfunction, xyz::AbstractVector{T})
    coordinate_transformation(u::AbstractSphericalEigenfunction, θ::T, ϕ::T)
    coordinate_transformation(u::AbstractSphericalEigenfunction, (θ::T, ϕ::T))
> Perform a coordinate transformation from the coordinate system used
  by the domain to that used by the eigenfunction.

  The point can be given in either Cartesian or spherical coordinates
  and the returned value will always be of the same type.
"""
function coordinate_transformation(u::AbstractSphericalEigenfunction,
                                   xyz::AbstractVector{T}
                                   ) where {T <: Union{arb, arb_series}}
    throw(ErrorException("coordinate_transform not implemented for eigenfunction of type $(typeof(u))"))
end

function coordinate_transformation(u::AbstractSphericalEigenfunction,
                                   θ::T,
                                   ϕ::T
                                   ) where {T <: Union{arb, arb_series}}
    spherical(coordinate_transformation(u, cartesian(θ, ϕ)))
end

function coordinate_transformation(u::AbstractSphericalEigenfunction,
                                   (θ, ϕ)::Union{Tuple{T, T},
                                                 NamedTuple{(:θ, :ϕ),Tuple{T, T}}},
                                   ) where {T <: Union{arb, arb_series}}
    coordinate_transformation(u, θ, ϕ)
end

"""
    u(xyz::AbstractVector{T}, λ::arb, k::Integer; boundary = nothing, notransform::Bool = false)
    u(θ::T, ϕ::T, λ::arb, k::Integer; boundary = nothing, notransform::Bool = false)
    u((θ::T, ϕ::T), λ::arb, k::Integer; boundary = nothing, notransform::Bool = false)
> Evaluate the k-th basis function for the eigenfunction

  The point can be given in either Cartesian or spherical coordinates.
  By default a coordinate transformation is applied to switch from the
  coordinate system used by the domain of u and that used by u. If
  notransform is true then do not perform this coordinate transform on
  the point first, this assumes that they already given in the
  coordinate system used by u.

  If boundary is set to an integer then the point is assumed to come
  from the corresponding boundary of the domain, some eigenfunctions
  are then identically equal to zero.

  See also: [`coordinate_transform`](@ref)
"""
function (u::AbstractSphericalEigenfunction)(xyz::AbstractVector{T},
                                             λ::arb,
                                             k::Integer;
                                             boundary = nothing,
                                             notransform::Bool = false
                                             ) where {T <: Union{arb, arb_series}}
    throw(ErrorException("evaluation of basis function not implemented"
                         *" for eigenfunction of type $(typeof(u))"))
end

function (u::AbstractSphericalEigenfunction)(θ::T,
                                             ϕ::T,
                                             λ::arb,
                                             k::Integer;
                                             boundary = nothing,
                                             notransform::Bool = false
                                             ) where {T <: Union{arb, arb_series}}
    throw(ErrorException("evaluation of basis function not implemented"
                         *" for eigenfunction of type $(typeof(u))"))
end

function (u::AbstractSphericalEigenfunction)((θ, ϕ)::Union{Tuple{T, T},
                                                           NamedTuple{(:θ, :ϕ),Tuple{T, T}}},
                                             λ::arb,
                                             k::Integer;
                                             boundary = nothing,
                                             notransform::Bool = false
                                             ) where {T <: Union{arb, arb_series}}
    u(θ, ϕ, λ, k, boundary = boundary, notransform = notransform)
end

"""
    u(xyz::AbstractVector{T}, λ::arb; boundary = nothing, notransform::Bool = false)
    u(θ::T, ϕ::T, λ::arb; boundary = nothing, notransform::Bool = false)
    u((θ::T, ϕ::T), λ::arb; boundary = nothing, notransform::Bool = false)
> Evaluate the eigenfunction.

  The point can be given in either Cartesian or spherical coordinates.
  By default a coordinate transformation is applied to switch from the
  coordinate system used by the domain of u and that used by u. If
  notransform is true then do not perform this coordinate transform on
  the point first, this assumes that they already given in the
  coordinate system used by u.

  If boundary is set to an integer then the point is assumed to come
  from the corresponding boundary of the domain, some eigenfunctions
  are then identically equal to zero.

  See also: [`coordinate_transform`](@ref)
"""
function (u::AbstractSphericalEigenfunction)(xyz::AbstractVector{T},
                                             λ::arb;
                                             boundary = nothing,
                                             notransform::Bool = false
                                             ) where {T <: Union{arb, arb_series}}
    res = u.domain.parent(0)

    if !notransform
        xyz = coordinate_transformation(u, xyz)
    end

    for k in 1:length(u.coefficients)
        res += u.coefficients[k]*u(xyz, λ, k, boundary = boundary, notransform = true)
        if (T == arb && !isfinite(res)) || (T == arb_series && !isfinite(res[end]))
            return res
        end
    end

    res
end

function (u::AbstractSphericalEigenfunction)(θ::T,
                                             ϕ::T,
                                             λ::arb;
                                             boundary = nothing,
                                             notransform::Bool = false
                                             ) where {T <: Union{arb, arb_series}}
    res = u.domain.parent(0)

    if !notransform
        θ, ϕ = coordinate_transformation(u, θ, ϕ)
    end

    for k in 1:length(u.coefficients)
        res += u.coefficients[k]*u(θ, ϕ, λ, k, boundary = boundary, notransform = true)
        if (T == arb && !isfinite(res)) || (T == arb_series && !isfinite(res[end]))
            return res
        end
    end

    res
end

function (u::AbstractSphericalEigenfunction)((θ, ϕ)::Union{Tuple{T, T},
                                                           NamedTuple{(:θ, :ϕ),Tuple{T, T}}},
                                             λ::arb;
                                             boundary = nothing,
                                             notransform::Bool = false
                                             ) where {T <: Union{arb, arb_series}}
    u(θ, ϕ, λ, boundary = boundary, notransform = notransform)
end

"""
    norm2(u::AbstractSphericalEigenfunction,
          λ::arb,
          (a, b, c))
> Lower bound the norm squared of the eigenfunction on the spherical
  triangle given by the three vertices `a`, `b` and `c`.

  The bound is computed by lower bounding `u^2` on the boundary and
  using a version of the minimum principle. To be able to use the
  minimum principle the area of the triangle needs to be small enough
  compared to the eigenvalue.
"""
function norm2(u::AbstractSphericalEigenfunction,
               λ::arb,
               (a, b, c))
    # Check that the area is small enough
    area = sum(anglesfromvertices(a, b, c)) - u.domain.parent(π)
    # FIXME: Use the correct maximum area here
    maximumarea = 4*u.domain.parent(π)^2

    if !(area < maximumarea)
        return u.domain.parent(0)
    end

    m = u.domain.parent(Inf)
    for (v, w) in [(a, b), (b, c), (c, a)]
        f = t -> -u(normalize(v .+ t.*(w - v)), λ)^2

        M = -enclosemaximum(f,
                            u.domain.parent(0),
                            u.domain.parent(1),
                            evaltype = :taylor,
                            n = div(length(coefficients(u)), 4),
                            atol = 1e-10,
                            rtol = 1e-3,
                            maxevals = 1000)

        m = min(m, M)
        if !isfinite(m) || !(m > 0)
            break
        end
    end

    res = m*area

    if isfinite(res) && res > 0
        return res
    else
        return u.domain.parent(0)
    end
end

function norm(u::AbstractSphericalEigenfunction,
              λ::arb)
    a, b, c = subtriangle(u.domain, ratio = 0.5)
    triangles = partitiontriangle(a, b, c, iterations = 1)

    res = sum(norm2(u, λ, triangle) for triangle in triangles)

    sqrt(res)
end
