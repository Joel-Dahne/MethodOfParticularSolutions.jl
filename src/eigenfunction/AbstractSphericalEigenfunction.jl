"""
    coefficients(u::AbstractSphericalEigenfunction)
> Return the coefficients of the eigenfunction.
"""
function coefficients(u::AbstractSphericalEigenfunction)
    u.coefficients
end

"""
    set_eigenfunction!(u::AbstractSphericalEigenfunction, coefficients)
> Set the coefficients for the expansion of the eigenfunction.
"""
function set_eigenfunction!(u::AbstractSphericalEigenfunction,
                            coefficients::Vector)
    resize!(u.coefficients, length(coefficients))
    copy!(u.coefficients, u.domain.parent.(coefficients))
end

"""
    coordinate_transformation(u::AbstractSphericalEigenfunction)
    coordinate_transformation(u::AbstractSphericalEigenfunction, θ::arb, ϕ::arb)
> Return a coordinate transformation T for switching to the coordinate
  system assumed by `u`.

  If θ and ϕ are given perform the coordinate transformation directly
  on θ and ϕ and return the result. Depending on the eigenfunction
  this can be more efficient if the transformation is not used
  multiple times.
"""
function coordinate_transformation(u::AbstractSphericalEigenfunction)
    throw(ErrorException("coordinate_transform not implemented for eigenfunction of type $(typeof(u))"))
end

function coordinate_transformation(u::AbstractSphericalEigenfunction,
                                   θ::T,
                                   ϕ::T) where {T <: Union{arb, arb_series}}
    coordinate_transformation(u)(θ, ϕ)
end

function coordinate_transformation(u::AbstractSphericalEigenfunction,
                                   (θ, ϕ)::Union{Tuple{T, T},
                                                   NamedTuple{(:θ, :ϕ),Tuple{T, T}}}
                                   ) where {T <: Union{arb, arb_series}}
    coordinate_transformation(u)(θ, ϕ)
end

"""
    u(θ::arb, ϕ::arb, λ::arb, k::Integer; notransform::Bool = false)
    u((θ::arb, ϕ::arb), λ::arb, k::Integer; notransform::Bool = false)
> Evaluate the k-th basis function for the eigenfunction with the
  given λ on the point given by (θ, ϕ).

  If notransform is true then do not perform a coordinate transform on
  θ and ϕ first, this assumes that they already given in the
  coordinate system used by u.

  See also: [`coordinate_transform`](@ref)
"""
function (u::AbstractSphericalEigenfunction)(θ::T,
                                             ϕ::T,
                                             λ::arb,
                                             k::Integer;
                                             notransform::Bool = false) where {T <: Union{arb, arb_series}}
    throw(ErrorException("evaluation of basis function not implemented for eigenfunction of type $(typeof(u))"))
end

function (u::AbstractSphericalEigenfunction)((θ, ϕ)::Union{Tuple{arb, arb},
                                                           NamedTuple{(:θ, :ϕ),Tuple{arb, arb}}},
                                             λ::arb,
                                             k::Integer;
                                             notransform::Bool = false)
    u(θ, ϕ, λ, k, notransform = notransform)
end

"""
    u(θ::arb, ϕ::arb, λ::arb; notransform::Bool = false)
    u((θ::arb, ϕ::arb), λ::arb; notransform::Bool = false)
> Evaluate the eigenfunction with the given λ on the point given by
  (θ, ϕ).

  If notransform is true then do not perform a coordinate transform on
  θ and ϕ first, this assumes that they already given in the
  coordinate system used by u.

  See also: [`coordinate_transform`](@ref)
"""
function (u::AbstractSphericalEigenfunction)(θ::T,
                                             ϕ::T,
                                             λ::arb;
                                             notransform::Bool = false) where {T <: Union{arb, arb_series}}
    res = u.domain.parent(0)

    if !notransform
        θ, ϕ = coordinate_transformation(u, θ, ϕ)
    end

    for k in 1:length(u.coefficients)
        res += u.coefficients[k]*u(θ, ϕ, λ, k, notransform = true)
    end

    res
end

function (u::AbstractSphericalEigenfunction)((θ, ϕ)::Union{Tuple{T, T},
                                                           NamedTuple{(:θ, :ϕ),Tuple{T, T}}},
                                             λ::arb;
                                             notransform::Bool = false) where {T <: Union{arb, arb_series}}
    u(θ, ϕ, λ, notransform = notransform)
end

function norm(u::AbstractSphericalEigenfunction,
              λ::arb)
    @error "no rigorous implementation of norm for $(typeof(u)), computing approximate norm"
    interior = interior_points(u.domain, 1000)
    sqrt(area(u.domain)*sum(abs(u(i, λ))^2 for i in interior)/length(interior))
end
