"""
    active_boundaries(u::AbstractSphericalEigenfunction)
> Return the boundaries which are active for the current
  eigenfunction. The eigenfunction is guaranteed to be identically
  equal to zero on the inactive boundaries.
"""
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
