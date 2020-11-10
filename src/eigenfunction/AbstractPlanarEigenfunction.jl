"""
    u(xy::AbstractVector{T}, λ::arb; boundary = nothing, notransform::Bool = false)
    u(r::T, θ::T, λ::arb; boundary = nothing, notransform::Bool = false)

Evaluate the eigenfunction.

The point can be given in either Cartesian or Polar coordinates. By
default a coordinate transformation is applied to switch from the
standard coordinate system to that used by `u`. If `notransform` is
true then this coordinate transformation is not performed, this
assumes that they already given in the coordinate system used by `u`.

If `boundary` is set to an integer then the point is assumed to come
from the corresponding boundary of the domain, some eigenfunctions are
then identically equal to zero.

See also: [`coordinate_transform`](@ref)
"""
function (u::AbstractPlanarEigenfunction) end

function (u::AbstractPlanarEigenfunction)(xy::AbstractVector{T},
                                          λ::arb;
                                          boundary = nothing,
                                          notransform::Bool = false
                                          ) where {T <: Union{arb, arb_series}}
    if !notransform
        xy = coordinate_transformation(u, xy)
    end

    coeffs = coefficients(u)
    res = zero(λ)
    for k in 1:length(coeffs)
        res += coeffs[k]*u(xy, λ, k, boundary = boundary, notransform = true)
        if (T == arb && !isfinite(res)) || (T == arb_series && !isfinite(res[end]))
            return res
        end
    end

    return res
end

function (u::AbstractPlanarEigenfunction)(r::T,
                                          θ::T,
                                          λ::arb;
                                          boundary = nothing,
                                          notransform::Bool = false
                                          ) where {T <: Union{arb, arb_series}}
    if !notransform
        r, θ = coordinate_transformation(u, r, θ)
    end

    coeffs = coefficients(u)
    res = zero(λ)
    for k in 1:length(coeffs)
        res += coeffs[k]*u(r, θ, λ, k, boundary = boundary, notransform = true)
        if (T == arb && !isfinite(res)) || (T == arb_series && !isfinite(res[end]))
            return res
        end
    end

    return res
end
