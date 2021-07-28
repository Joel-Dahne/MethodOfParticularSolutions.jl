function SphericalInteriorEigenfunction(
    domain::SphericalTriangle,
    θ::arb,
    ϕ::arb;
    stride::Int = 1,
)
    SphericalInteriorEigenfunction(domain, θ, ϕ, stride, arb[])
end

function SphericalInteriorEigenfunction(
    domain::SphericalTriangle,
    xyz::AbstractVector{arb};
    stride::Int = 1,
)
    SphericalInteriorEigenfunction(domain, spherical(xyz)..., stride = stride)
end

function SphericalInteriorEigenfunction(domain::SphericalTriangle; stride::Int = 1)
    SphericalInteriorEigenfunction(domain, spherical(center(domain))..., stride = stride)
end

function Base.show(io::IO, u::SphericalInteriorEigenfunction)
    println(io, "Interior eigenfunction")
    if !haskey(io, :compact) || !io[:compact]
        println(io, "interior point: (θ, ϕ) = ($(u.θ), $(u.ϕ))")
        println(io, "domain: $(u.domain)")
        print(io, "number of set coefficients: $(length(u.coefficients))")
    end
end

function recompute!(u::SphericalInteriorEigenfunction)
    u.θ, u.ϕ = spherical(center(u.domain))
    u
end

function coordinate_transformation(
    u::SphericalInteriorEigenfunction,
    xyz::AbstractVector{T},
) where {T<:Union{arb,arb_series}}
    # Rotate by u.ϕ along the along the z-axis so that the ϕ value of
    # the point becomes zero, the rotate by u.θ along the y-axis so
    # that the point ends up on the north pole.
    L = LinearMap(RotYZ(-u.θ, -u.ϕ))
    L(xyz)
end

function coordinate_transformation(
    u::SphericalInteriorEigenfunction,
    θ::T,
    ϕ::T,
) where {T<:Union{arb,arb_series}}
    # TODO: The performance could most likely be improved by
    # performing the rotation the z-axis while still in spherical
    # coordinates.
    spherical(coordinate_transformation(u, cartesian(θ, ϕ)))
end

function (u::SphericalInteriorEigenfunction)(
    xyz::AbstractVector{T},
    λ,
    k::Integer;
    boundary = nothing,
    notransform::Bool = false,
) where {T<:Union{arb,arb_series}}
    λ = u.parent(λ)

    k = 1 + (k - 1) * u.stride

    ν::arb = -0.5 + sqrt(0.25 + λ)
    μ::arb = u.domain.parent(div(k, 2))
    if !notransform
        xyz = coordinate_transformation(u, xyz)
    end

    if k == 1
        return legendre_p_safe(ν, μ, xyz[3])
    elseif k % 2 == 0
        ϕ = atan(xyz[2], xyz[1])
        return legendre_p_safe(ν, μ, xyz[3]) * sin(μ * ϕ)
    else
        ϕ = atan(xyz[2], xyz[1])
        return legendre_p_safe(ν, μ, xyz[3]) * cos(μ * ϕ)
    end
end


function (u::SphericalInteriorEigenfunction)(
    θ::T,
    ϕ::T,
    λ,
    k::Integer;
    boundary = nothing,
    notransform::Bool = false,
) where {T<:Union{arb,arb_series}}
    λ = u.parent(λ)

    k = 1 + (k - 1) * u.stride

    ν::arb = -0.5 + sqrt(0.25 + λ)
    μ::arb = u.domain.parent(div(k, 2))
    if !notransform
        θ, ϕ = coordinate_transformation(u, θ, ϕ)
    end
    if k == 1
        return legendre_p_safe(ν, μ, cos(θ))
    elseif k % 2 == 0
        return legendre_p_safe(ν, μ, cos(θ)) * sin(μ * ϕ)
    else
        return legendre_p_safe(ν, μ, cos(θ)) * cos(μ * ϕ)
    end
end
