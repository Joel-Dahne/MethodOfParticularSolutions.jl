LShapeEigenfunction(domain::LShape; stride::Int = 1) =
    LShapeEigenfunction(domain, stride, arb[])

function Base.show(io::IO, u::LShapeEigenfunction)
    println(io, "Eigenfunction for the L-shaped domain")
    if !haskey(io, :compact) || !io[:compact]
        println(io, "domain: $(u.domain)")
        print(io, "number of set coefficients: $(length(u.coefficients))")
    end
end

function Base.getproperty(u::LShapeEigenfunction, name::Symbol)
    if name == :parent
        return u.domain.parent
    else
        return getfield(u, name)
    end
end

active_boundaries(::LShape, ::LShapeEigenfunction) = 1:4

function (u::LShapeEigenfunction)(
    (r, θ)::Tuple{T,T},
    λ::arb,
    ks::UnitRange{Int};
    boundary = nothing,
) where {T<:Union{arb,arb_series}}
    rsqrtλ = r * sqrt(λ)
    res = similar(ks, T)
    for i in eachindex(ks)
        k = 1 + (ks[i] - 1) * u.stride

        ν = u.parent(2k // 3)

        res[i] = bessel_j(ν, rsqrtλ) * sin(ν * θ)
    end

    return res
end

function (u::LShapeEigenfunction)(
    (r, θ)::Tuple{T,T},
    λ::arb;
    boundary = nothing,
) where {T<:Union{arb,arb_series}}
    coeffs = coefficients(u)
    return sum(coeffs .* u((r, θ), λ, 1:length(coeffs), boundary = boundary))
end

function norm(domain::LShape, u::LShapeEigenfunction, λ::arb)
    θ_integral = 3 // 4 * domain.parent(π)

    CC = ComplexField(domain.parent.prec)
    # The integrals goes from zero to the lower bound for θ. However
    # the function has a branch cut at zero and Arb has problem
    # handling this. We therefore integrate a small distance away from
    # zero.
    a = CC(1e-1)
    b = CC(1)

    r_integral = domain.parent(0)
    for k = 1:min(4, length(u.coefficients))
        ν::arb = domain.parent(2k // 3)
        c2 = u.coefficients[k]^2

        f = r -> begin
            r = real(r)
            if contains_negative(r)
                return CC(NaN)
            end

            CC(c2 * r * bessel_j(ν, sqrt(λ) * r)^2)
        end
        y = Nemo.integrate(CC, f, a, b)
        r_integral += real(y)
    end

    return sqrt(θ_integral * r_integral)
end
