abstract type AbstractEigenfunction end

abstract type AbstractVertexEigenfunction <: AbstractEigenfunction end

mutable struct SphericalVertexEigenfunction{T <: Union{fmpq, arb}} <: AbstractVertexEigenfunction
    domain::SphericalTriangle{T}
    vertex::Int
    coefficients::Vector{arb}
end
