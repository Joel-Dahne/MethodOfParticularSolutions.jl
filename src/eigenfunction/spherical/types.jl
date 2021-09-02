abstract type AbstractSphericalEigenfunction{S,T} <: AbstractEigenfunction{S,T} end

mutable struct SphericalVertexEigenfunction{T<:Union{fmpq,arb}} <:
               AbstractSphericalEigenfunction{arb,T}
    domain::SphericalTriangle{T}
    vertex::Int
    stride::Int
    coefficients::Vector{arb}
end

mutable struct SphericalInteriorEigenfunction{T<:Union{fmpq,arb}} <:
               AbstractSphericalEigenfunction{arb,T}
    domain::SphericalTriangle{T}
    θ::arb # θ value for the interior point
    ϕ::arb # ϕ value for the interior point
    stride::Int
    coefficients::Vector{arb}
end

mutable struct SphericalCombinedEigenfunction{T<:Union{fmpq,arb}} <:
               AbstractSphericalEigenfunction{arb,T}
    domain::SphericalTriangle{T}
    us::Vector{<:AbstractSphericalEigenfunction}
    orders::Vector{Int}
end

mutable struct KrewerasEigenfunction{T<:Union{fmpq,arb}} <:
               AbstractSphericalEigenfunction{arb,T}
    domain::SphericalTriangle{T}
    coefficients::Vector{arb}
end