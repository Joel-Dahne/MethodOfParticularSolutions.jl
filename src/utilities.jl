function cartesian_from_polar(r, θ)
    s, c = sincos(θ)
    x = r*c
    y = r*s
    return SVector(x, y)
end

function polar_from_cartesian(xy)
    r = sqrt(xy[1]^2 + xy[2]^2)
    θ = atan(xy[2], xy[1])
    return r, θ
end

function cartesian(θ, ϕ)
    sθ, cθ = sincos(θ)
    sϕ, cϕ = sincos(ϕ)
    SVector(sθ*cϕ, sθ*sϕ, cθ)
end

cartesian((θ, ϕ)) = cartesian(θ, ϕ)

function spherical(xyz)
    (θ = acos(xyz[3]), ϕ = atan(xyz[2], xyz[1]))
end

function Base.Float64((θ, ϕ)::NamedTuple{(:θ, :ϕ),Tuple{arb,arb}})
    (θ = Float64(θ), ϕ = Float64(ϕ))
end

function LinearAlgebra.dot(xyz1::SVector{3, arb}, xyz2::SVector{3, arb})
    sum(xyz1.*xyz2)
end

LinearAlgebra.norm(x::SVector{N,arb}) where {N} = sqrt(sum(x.^2))

LinearAlgebra.normalize(x::SVector{N,arb}) where {N} = x./sqrt(sum(x.^2))

LinearAlgebra.normalize(x::SVector{N,arb_series}) where {N} = x./sqrt(sum(x.^2))
