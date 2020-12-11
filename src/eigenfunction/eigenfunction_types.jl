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

# TODO: Currently this can't handle the case when some eigenfunctions
# are identically equal to zero on some parts of the boundary. As long
# as we are working with lightning and interior eigenfunctions this is
# fine, but it means it's not very useful for vertex eigenfunctions.
struct LinkedEigenfunction{T<:AbstractPlanarEigenfunction} <: AbstractPlanarEigenfunction
    us::Vector{T}
    extra_coefficients::Vector{arb}

    function LinkedEigenfunction(
        us::Vector{T},
        extra_coefficients::Vector = ones(length(us)),
    ) where {T <: AbstractPlanarEigenfunction}
        isempty(us) && throw(ArgumentError("us must not be empty"))
        length(us) == length(extra_coefficients) ||
            throw(ArgumentError("us and extra_coefficients must have the same size"))

        extra_coefficients = first(us).parent.(extra_coefficients)

        return new{T}(us, extra_coefficients)
    end
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
