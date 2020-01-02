abstract type AbstractEigenfunction end

abstract type AbstractVertexEigenfunction <: AbstractEigenfunction end

struct SphericalVertexEigenfunction{T <: Union{fmpq, arb}} <: AbstractVertexEigenfunction
    domain::SphericalTriangle{T}
    vertex::Int
    coefficients::Vector{arb}
end
