function SphericalInteriorEigenfunction(domain::SphericalTriangle,
                                        θ::arb,
                                        ϕ::arb)
    SphericalInteriorEigenfunction(domain, θ, ϕ, arb[])
end

function Base.show(io::IO, u::SphericalInteriorEigenfunction)
    println(io, "Interior eigenfunction")
    if !haskey(io, :compact) || !io[:compact]
        println(io, "interior point: (θ, ϕ) = ($(u.θ), $(u.ϕ))")
        println(io, "domain: $(u.domain)")
        print(io, "number of set coefficients: $(length(u.coefficients))")
    end
end

"""
    coordinate_transformation(u::SphericalInteriorEigenfunction)
> Return a coordinate transformation T which switches from spherical
  coordinates (θ, ϕ) to (θ', ϕ'), in the new coordinate system the
  point from which the eigenfunction originates is put at the north
  pole.

  This is computed by switching to cartesian coordinates, performing
  two rotations and then switching back to spherical coordinates.
"""
function coordinate_transformation(u::SphericalInteriorEigenfunction)
    # TODO: The performance could most likely be improved

    # Rotate by u.ϕ along the along the z-axis so that the ϕ value of
    # the point becomes zero, the rotate by u.θ along the y-axis so
    # that the point ends up on the north pole.
    L = LinearMap(RotYZ(-u.θ, -u.ϕ))
    T = (θ, ϕ) -> begin
        xyz = cartesian(θ, ϕ)
        xyz = L(xyz)
        (θ, ϕ) = spherical(xyz)
        return θ, ϕ
    end

    return T
end

function (u::SphericalInteriorEigenfunction)(θ::T,
                                             ϕ::T,
                                             λ::arb,
                                             k::Integer;
                                             notransform::Bool = false) where {T <: Union{arb, arb_series}}
    ν = u.domain.parent(-0.5) + sqrt(u.domain.parent(0.25) + λ)
    if !notransform
        θ, ϕ = coordinate_transformation(u, θ, ϕ)
    end
    if k == 1
        μ = u.domain.parent(0)
        return legendre_p_safe(ν, μ, cos(θ))
    elseif k % 2 == 0
        μ = u.domain.parent(div(k, 2))
        return legendre_p_safe(ν, μ, cos(θ))*sin(μ*ϕ)
    else
        μ = u.domain.parent(div(k, 2))
        return legendre_p_safe(ν, μ, cos(θ))*cos(μ*ϕ)
    end
end
