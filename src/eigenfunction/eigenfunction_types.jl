abstract type AbstractEigenfunction end

abstract type AbstractSphericalEigenfunction <: AbstractEigenfunction end

struct SphericalVertexEigenfunction{T <: Union{fmpq, arb}} <: AbstractSphericalEigenfunction
    domain::SphericalTriangle{T}
    vertex::Int
    coefficients::Vector{arb}
end

struct SphericalInteriorEigenfunction{T <: Union{fmpq, arb}} <: AbstractSphericalEigenfunction
    domain::SphericalTriangle{T}
    θ::arb # θ value for the interior point
    ϕ::arb # ϕ value for the interior point
    coefficients::Vector{arb}
end
