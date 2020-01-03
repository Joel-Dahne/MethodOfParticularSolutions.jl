function cartesian(θ, ϕ)
    [sin(θ)*cos(ϕ), sin(θ)*sin(ϕ), cos(θ)]
end

cartesian((θ, ϕ)) = cartesian(θ, ϕ)

function spherical(xyz)
    r = sqrt(sum(xyz.^2))
    (θ = acos(xyz[3]/r), ϕ = atan(xyz[2], xyz[1]))
end

function Base.Float64((θ, ϕ)::NamedTuple{(:θ, :ϕ),Tuple{arb,arb}})
    (θ = Float64(θ), ϕ = Float64(ϕ))
end
