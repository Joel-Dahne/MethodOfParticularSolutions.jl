function LShapeEigenfunction(domain::LShape;
                             stride::Int = 1)
    LShapeEigenfunction(domain, stride, arb[])
end

function Base.show(io::IO, u::LShapeEigenfunction)
    println(io, "Eigenfunction for the L-shaped domain")
    if !haskey(io, :compact) || !io[:compact]
        println(io, "domain: $(u.domain)")
        print(io, "number of set coefficients: $(length(u.coefficients))")
    end
end

function (u::LShapeEigenfunction)(r::T,
                                  θ::T,
                                  λ::arb,
                                  k::Integer;
                                  boundary = nothing
                                  ) where {T <: Union{arb, arb_series}}
    k = 1 + (k - 1)*u.stride

    ν = u.domain.parent(2k//3)

    bessel_j(ν, sqrt(λ)*r)*sin(ν*θ)
end

function (u::LShapeEigenfunction)((r, θ)::Tuple{T, T},
                                  λ::arb,
                                  k::Integer;
                                  boundary = nothing
                                  ) where {T <: Union{arb, arb_series}}
    u(r, θ, λ, k, boundary = boundary)
end

function (u::LShapeEigenfunction)(r::T,
                                  θ::T,
                                  λ::arb;
                                  boundary = nothing
                                  ) where {T <: Union{arb, arb_series}}
    res = u.domain.parent(0)

    for k in 1:length(u.coefficients)
        res += u.coefficients[k]*u(r, θ, λ, k, boundary = boundary)
    end

    res
end

function (u::LShapeEigenfunction)((r, θ)::Tuple{T, T},
                                  λ::arb;
                                  boundary = nothing
                                  ) where {T <: Union{arb, arb_series}}
    u(r, θ, λ, boundary = boundary)
end

function norm(u::LShapeEigenfunction,
              λ::arb)
    θ_integral = 3//4*u.domain.parent(π)

    CC = ComplexField(u.domain.parent.prec)
    # The integrals goes from zero to the lower bound for θ. However
    # the function has a branch cut at zero and Arb has problem
    # handling this. We therefore integrate a small distance away from
    # zero.
    a = CC(1e-1)
    b = CC(1)

    r_integral = u.domain.parent(0)
    for k in 1:min(4, length(u.coefficients))
        ν::arb = u.domain.parent(2k//3)
        c2 = u.coefficients[k]^2

        f = r -> begin
            r = real(r)
            if contains_negative(r)
                return CC(NaN)
            end

            CC(c2*r*bessel_j(ν, sqrt(λ)*r)^2)
        end
        y = Nemo.integrate(CC, f, a, b)
        r_integral += real(y)
    end

    sqrt(θ_integral*r_integral)
end
