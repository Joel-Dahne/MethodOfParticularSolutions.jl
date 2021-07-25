"""
    u(xy::AbstractVector, λ; boundary = nothing, notransform::Bool = false)
    u(xy::AbstractVector, λ, k::Integer; boundary = nothing, notransform::Bool = false)
    u(xy::AbstractVector, λ, ks::UnitRange{Int}; boundary = nothing, notransform::Bool = false)

Evaluate the eigenfunction at the point `xy` given in Cartesian
coordinates.

By default a coordinate transformation is applied to switch from the
standard coordinate system to that used by `u`. If `notransform` is
true then this coordinate transformation is not performed, this
assumes that they already given in the coordinate system used by `u`.

If `k::Integer` is given then evaluate only that basis function and
don't multiply with any coefficient. If `ks::UnitRange{Int}` is given
then return a vector with the values of the corresponding basis
functions.

If `boundary` is set to an integer then the point is assumed to come
from the corresponding boundary of the domain, some eigenfunctions are
then identically equal to zero.

See also: [`coordinate_transform`](@ref)
"""
function (u::AbstractPlanarEigenfunction) end

(u::AbstractPlanarEigenfunction)(
    xy::AbstractVector,
    λ,
    k::Integer;
    boundary = nothing,
    notransform::Bool = false,
) = u(xy, λ, k:k; boundary, notransform)

function (u::AbstractPlanarEigenfunction{S,T})(
    xy::AbstractVector,
    λ;
    boundary = nothing,
    notransform::Bool = false,
) where {S,T}
    coeffs = coefficients(u)

    if isempty(coeffs)
        # If there are no coefficients we can't rely on the evaluation
        # of the eigenfunction to determine the return type.
        if S == arb
            if eltype(xy) == arb_series
                return arb_series(parent(first(xy).poly)(0), length(first(xy)))
            else
                return u.parent(0)
            end
        else
            return zero(S)
        end
    end

    return sum(coeffs .* u(xy, λ, 1:length(coeffs); boundary, notransform))
end
