function cartesian_from_polar(r, θ)
    s, c = sincos(θ)
    x = r * c
    y = r * s
    return SVector(x, y)
end

# TODO: tighter enclosures for wide intervals
polar_from_cartesian(xy::AbstractVector) = (hypot(xy[1], xy[2]), atan(xy[2], xy[1]))

function polar_from_cartesian(xy::AbstractVector{arb_series})
    r = sqrt(xy[1]^2 + xy[2]^2)

    if contains_zero(xy[2][0]) && !iszero(xy[2][0]) && !isnonnegative(xy[1][0])
        θ = atan(xy[2], -xy[1]) + base_ring(xy[1].poly)(π)
    else
        θ = atan(xy[2], xy[1])
    end

    return r, θ
end

function polar_from_cartesian(xy::AbstractVector{arb})
    r = hypot(xy[1], xy[2])

    if contains_zero(xy[2]) && !iszero(xy[2]) && !isnonnegative(xy[1])
        θ = atan(xy[2], -xy[1]) + parent(xy[1])(π)
    else
        θ = atan(xy[2], xy[1])
    end

    return r, θ
end

function cartesian(θ, ϕ)
    sθ, cθ = sincos(θ)
    sϕ, cϕ = sincos(ϕ)
    SVector(sθ * cϕ, sθ * sϕ, cθ)
end

cartesian((θ, ϕ)) = cartesian(θ, ϕ)

function spherical(xyz)
    (θ = acos(xyz[3]), ϕ = atan(xyz[2], xyz[1]))
end

function Base.Float64((θ, ϕ)::NamedTuple{(:θ, :ϕ),Tuple{arb,arb}})
    (θ = Float64(θ), ϕ = Float64(ϕ))
end

function LinearAlgebra.dot(xyz1::SVector{N,arb}, xyz2::SVector{N,arb}) where {N}
    sum(xyz1 .* xyz2)
end

LinearAlgebra.norm(x::SVector{N,arb}) where {N} = sqrt(sum(x .^ 2))

LinearAlgebra.normalize(x::SVector{N,arb}) where {N} = x ./ sqrt(sum(x .^ 2))

LinearAlgebra.normalize(x::SVector{N,arb_series}) where {N} = x ./ sqrt(sum(x .^ 2))
