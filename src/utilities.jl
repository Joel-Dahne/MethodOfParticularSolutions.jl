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

function LinearAlgebra.norm(xyz::SVector{3, arb})
    sqrt(sum(xyz.^2))
end

function LinearAlgebra.normalize(xyz::SVector{3, arb})
    r = sqrt(sum(xyz.^2))
    xyz./r
end

function LinearAlgebra.normalize(xyz::SVector{3, arb_series})
    r = sqrt(sum(xyz.^2))
    xyz./r
end
