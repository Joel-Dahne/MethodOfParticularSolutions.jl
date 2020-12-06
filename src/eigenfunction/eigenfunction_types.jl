abstract type AbstractEigenfunction end

abstract type AbstractSphericalEigenfunction <: AbstractEigenfunction end

abstract type AbstractPlanarEigenfunction <: AbstractEigenfunction end

###
### Spherical eigenfunctions
###

mutable struct SphericalVertexEigenfunction{T <: Union{fmpq, arb}} <: AbstractSphericalEigenfunction
    domain::SphericalTriangle{T}
    vertex::Int
    stride::Int
    coefficients::Vector{arb}
end

mutable struct SphericalInteriorEigenfunction{T <: Union{fmpq, arb}} <: AbstractSphericalEigenfunction
    domain::SphericalTriangle{T}
    θ::arb # θ value for the interior point
    ϕ::arb # ϕ value for the interior point
    stride::Int
    coefficients::Vector{arb}
end

mutable struct SphericalCombinedEigenfunction{T <: Union{fmpq, arb}} <: AbstractSphericalEigenfunction
    domain::SphericalTriangle{T}
    us::Vector{<:AbstractSphericalEigenfunction}
    orders::Vector{Int}
end

mutable struct KrewerasEigenfunction{T <: Union{fmpq, arb}} <: AbstractSphericalEigenfunction
    domain::SphericalTriangle{T}
    coefficients::Vector{arb}
end

###
### Planar eigenfunctions
###

struct StandaloneVertexEigenfunction{T <: Union{fmpq, arb}} <: AbstractPlanarEigenfunction
    vertex::SVector{2,arb}
    orientation::T
    θ::T
    stride::Int
    coefficients::Vector{arb}
    parent::ArbField
end

struct StandaloneInteriorEigenfunction <: AbstractPlanarEigenfunction
    vertex::SVector{2,arb}
    stride::Int
    coefficients::Vector{arb}
    parent::ArbField
end

struct StandaloneLightningEigenfunction{T <: Union{fmpq, arb}} <: AbstractPlanarEigenfunction
    vertex::SVector{2,arb}
    orientation::T
    θ::T
    l::arb
    σ::arb
    coefficients::Vector{arb}
    parent::ArbField
end

mutable struct CombinedEigenfunction <: AbstractPlanarEigenfunction
    domain::AbstractPlanarDomain
    us::Vector{<:AbstractPlanarEigenfunction}
    boundary_to_us::OrderedDict{Int,BitSet}
    orders::Vector{Int}
end

mutable struct LShapeEigenfunction <: AbstractPlanarEigenfunction
    domain::LShape
    stride::Int
    coefficients::Vector{arb}
end
