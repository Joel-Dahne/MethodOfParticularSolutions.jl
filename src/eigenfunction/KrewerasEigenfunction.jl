function KrewerasEigenfunction(domain::SphericalTriangle)
    if domain.angles[1] != fmpq(2//3) || domain.angles[2] != fmpq(2//3) || domain.angles[3] != fmpq(2//3)
        throw(ArgumentError("domain must be the Kreweras triangle with angles (2π/3, 2π/3, 2π/3)"))
    end
    KrewerasEigenfunction(domain, arb[])
end

function Base.show(io::IO, u::KrewerasEigenfunction)
    println(io, "Kreweras eigenfunction")
    if !haskey(io, :compact) || !io[:compact]
        println(io, "domain: $(u.domain)")
        print(io, "number of set coefficients: $(length(u.coefficients))")
    end
end

function active_boundaries(u::KrewerasEigenfunction)
    (true, false, false)
end

function (u::KrewerasEigenfunction)(xyz::AbstractVector{T},
                                    λ::arb,
                                    k::Integer;
                                    boundary = nothing,
                                    notransform::Bool = false
                                    ) where {T <: Union{arb, arb_series}}
    ν::arb = -0.5 + sqrt(0.25 + λ)
    if k%2 == 1
        μ::arb = u.domain.parent(-k*inv(u.domain.angles[1]))

        # Create a linear map for rotating so that the next vertex ends up
        # at the north pole
        α = -acos(vertex(u.domain, 2)[3])
        β = angle(u.domain, 2) - u.domain.parent(π)
        L = LinearMap(RotZY(β, α))

        res = u.domain.parent(0)

        if isnothing(boundary) || boundary == 1
            ϕ = atan(xyz[2], xyz[1])
            res += legendre_p_safe(ν, μ, xyz[3])*sin(μ*ϕ)
        end

        xyz = L(xyz)
        if isnothing(boundary) || boundary == 2
            ϕ = atan(xyz[2], xyz[1])
            res += legendre_p_safe(ν, μ, xyz[3])*sin(μ*ϕ)
        end

        if isnothing(boundary) || boundary == 3
            xyz = L(xyz)
            ϕ = atan(xyz[2], xyz[1])
            res += legendre_p_safe(ν, μ, xyz[3])*sin(μ*ϕ)
        end

        return res
    else
        μ = u.domain.parent(3div(k - 2, 2))

        (θ, φ) = spherical(center(u.domain))
        L = LinearMap(RotYZ(-θ, -φ))
        xyz = L(xyz)
        ϕ = atan(xyz[2], xyz[1])
        return legendre_p_safe(ν, μ, xyz[3])*cos(μ*ϕ)
    end
end

function (u::KrewerasEigenfunction)(xyz::AbstractVector{T},
                                    λ::arb;
                                    boundary = nothing,
                                    notransform::Bool = false
                                    ) where {T <: Union{arb, arb_series}}
    res = u.domain.parent(0)

    for k in 1:length(u.coefficients)
        res += u.coefficients[k]*u(xyz, λ, k, boundary = boundary)
        if (T == arb && !isfinite(res)) || (T == arb_series && !isfinite(res[end]))
            return res
        end
    end

    res
end
