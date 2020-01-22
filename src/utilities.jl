function cartesian(θ, ϕ)
    sθ, cθ = sincos(θ)
    sϕ, cϕ = sincos(ϕ)
    [sθ*cϕ, sθ*sϕ, cθ]
end

function cartesian(θ::arb_poly, ϕ::arb_poly)
    n = length(θ)
    sθ, cθ = sincos(θ)
    sϕ, cϕ = sincos(ϕ)
    [mullow(sθ, cϕ, n), mullow(sθ, sϕ, n), cθ]
end

cartesian((θ, ϕ)) = cartesian(θ, ϕ)

function spherical(xyz)
    r = sqrt(sum(xyz.^2))
    (θ = acos(xyz[3]/r), ϕ = atan(xyz[2], xyz[1]))
end

function Base.Float64((θ, ϕ)::NamedTuple{(:θ, :ϕ),Tuple{arb,arb}})
    (θ = Float64(θ), ϕ = Float64(ϕ))
end
