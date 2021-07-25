"""
    u(xy::AbstractVector{T}, λ::arb, k::Int; boundary = nothing, notransform::Bool = false)
    u(xy::AbstractVector{T}, λ::arb, ks::UnitRange{Int}; boundary = nothing, notransform::Bool = false)
    u(xy::AbstractVector{T}, λ::arb; boundary = nothing, notransform::Bool = false)

Evaluate the eigenfunction.

The point can be given in either Cartesian or Polar coordinates. By
default a coordinate transformation is applied to switch from the
standard coordinate system to that used by `u`. If `notransform` is
true then this coordinate transformation is not performed, this
assumes that they already given in the coordinate system used by `u`.

If `k::Int` is given then evaluate only that basis function and don't
multiply with any coefficient. If `ks::UnitRange{Int}` is given then
return a vector with the values of the corresponding basis functions.

If `boundary` is set to an integer then the point is assumed to come
from the corresponding boundary of the domain, some eigenfunctions are
then identically equal to zero.

See also: [`coordinate_transform`](@ref)
"""
function (u::AbstractPlanarEigenfunction) end

(u::AbstractPlanarEigenfunction)(
    xy::AbstractVector{T},
    λ::arb,
    k::Integer;
    boundary = nothing,
    notransform::Bool = false,
) where {T<:Union{arb,arb_series}} =
    u(xy, λ, k:k, boundary = boundary, notransform = notransform)

function (u::AbstractPlanarEigenfunction)(
    xy::AbstractVector{T},
    λ::arb;
    boundary = nothing,
    notransform::Bool = false,
) where {T<:Union{arb,arb_series}}
    # TODO: This computes all terms even if some terms become NaN very
    # early on. Could potentially be optimized.
    coeffs = coefficients(u)
    isempty(coeffs) && return zero(λ)
    return sum(coeffs .* u(xy, λ, 1:length(coeffs), boundary = boundary))
end
