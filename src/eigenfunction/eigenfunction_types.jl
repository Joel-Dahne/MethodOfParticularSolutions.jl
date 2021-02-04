abstract type AbstractEigenfunction end

abstract type AbstractSphericalEigenfunction <: AbstractEigenfunction end

abstract type AbstractPlanarEigenfunction <: AbstractEigenfunction end

###
### Spherical eigenfunctions
###

mutable struct SphericalVertexEigenfunction{T<:Union{fmpq,arb}} <:
               AbstractSphericalEigenfunction
    domain::SphericalTriangle{T}
    vertex::Int
    stride::Int
    coefficients::Vector{arb}
end

mutable struct SphericalInteriorEigenfunction{T<:Union{fmpq,arb}} <:
               AbstractSphericalEigenfunction
    domain::SphericalTriangle{T}
    θ::arb # θ value for the interior point
    ϕ::arb # ϕ value for the interior point
    stride::Int
    coefficients::Vector{arb}
end

mutable struct SphericalCombinedEigenfunction{T<:Union{fmpq,arb}} <:
               AbstractSphericalEigenfunction
    domain::SphericalTriangle{T}
    us::Vector{<:AbstractSphericalEigenfunction}
    orders::Vector{Int}
end

mutable struct KrewerasEigenfunction{T<:Union{fmpq,arb}} <: AbstractSphericalEigenfunction
    domain::SphericalTriangle{T}
    coefficients::Vector{arb}
end

###
### Planar eigenfunctions
###

struct StandaloneVertexEigenfunction{T<:Union{fmpq,arb}} <: AbstractPlanarEigenfunction
    vertex::SVector{2,arb}
    orientation::T
    θ::T
    stride::Int
    reversed::Bool
    coefficients::Vector{arb}
    parent::ArbField

    function StandaloneVertexEigenfunction(
        vertex::SVector{2,arb},
        orientation::T,
        θ::T;
        stride::Integer = 1,
        reversed::Bool = false,
        parent::ArbField = parent(vertex[1]),
    ) where {T<:Union{Rational,arb,fmpq}}
        if T == fmpq || T <: Rational
            orientation = fmpq(orientation)
            S = fmpq
        else
            S = arb
        end
        return new{S}(vertex, orientation, θ, stride, reversed, arb[], parent)
    end
end

""""
    StandaloneInteriorEigenfunction(vertex, orientation, parent; stride = 1, even = false)

An eigenfunction consisting of the functions `bessel_j(ν,
r*√λ)*sin(ν*θ)` and `bessel_j(ν, ν*√λ)*cos(j*θ)` for `ν` = 0, 1, 2,
... Here `r` and `θ` are polar coordinates around `vertex` take so
that `θ = 0` at the given orientation.

The values of `ν` that are used are 0, `stride`, `2stride`,
`3stride`,...

If `even` is true then use only the function with `cos(j*θ)`.
"""
struct StandaloneInteriorEigenfunction{T<:Union{fmpq,arb}} <: AbstractPlanarEigenfunction
    vertex::SVector{2,arb}
    orientation::T
    stride::Int
    even::Bool
    coefficients::Vector{arb}
    parent::ArbField
end

struct StandaloneLightningEigenfunction{T<:Union{AbstractFloat,arb},S<:Union{fmpq,arb}} <:
       AbstractPlanarEigenfunction
    vertex::SVector{2,T}
    orientation::S
    θ::S
    l::T
    σ::T
    even::Bool
    reversed::Bool
    coefficients::Vector{T}
    parent::ArbField

    function StandaloneLightningEigenfunction(
        vertex::SVector{2,arb},
        orientation::S,
        θ::S,
        parent::ArbField = parent(vertex[1]);
        l::arb = parent(1),
        σ::arb = parent(4),
        even::Bool = false,
        reversed::Bool = false,
    ) where {S<:Union{arb,fmpq}}
        return new{arb,S}(vertex, orientation, θ, l, σ, even, reversed, arb[], parent)
    end

    function StandaloneLightningEigenfunction{arb,S}(
        vertex::SVector{2,arb},
        orientation::S,
        θ::S,
        parent::ArbField = parent(vertex[1]);
        l = 1,
        σ = 4,
        even::Bool = false,
        reversed::Bool = false,
    ) where {S<:Union{arb,fmpq}}
        return new{arb,S}(
            vertex,
            orientation,
            θ,
            parent(l),
            parent(σ),
            even,
            reversed,
            arb[],
            parent,
        )
    end

    function StandaloneLightningEigenfunction{T,S}(
        vertex::AbstractVector,
        orientation::S,
        θ::S;
        l = one(T),
        σ = 4one(T),
        even::Bool = false,
        reversed::Bool = false,
    ) where {T<:AbstractFloat,S<:Union{arb,fmpq}}
        return new{T,S}(
            convert(SVector{2,T}, vertex),
            orientation,
            θ,
            convert(T, l),
            convert(T, σ),
            even,
            reversed,
            T[],
            ArbField(precision(T)),
        )
    end
end

struct LinkedEigenfunction{T<:AbstractPlanarEigenfunction} <: AbstractPlanarEigenfunction
    us::Vector{T}
    extra_coefficients::Vector{arb}
    excluded_boundaries::Vector{BitSet}

    function LinkedEigenfunction(
        us::Vector{T},
        extra_coefficients::Vector = ones(length(us));
        excluded_boundaries::Vector = fill(BitSet(), length(us)),
    ) where {T<:AbstractPlanarEigenfunction}
        isempty(us) && throw(ArgumentError("us must not be empty"))
        length(us) == length(extra_coefficients) ||
            throw(ArgumentError("us and extra_coefficients must have the same size"))

        extra_coefficients = first(us).parent.(extra_coefficients)

        return new{T}(us, extra_coefficients, excluded_boundaries)
    end
end

mutable struct CombinedEigenfunction <: AbstractPlanarEigenfunction
    domain::AbstractPlanarDomain
    us::Vector{<:AbstractPlanarEigenfunction}
    boundary_to_us::OrderedDict{Int,BitSet}
    even_boundaries::BitSet
    orders::Vector{Int}
end

mutable struct LShapeEigenfunction <: AbstractPlanarEigenfunction
    domain::LShape
    stride::Int
    coefficients::Vector{arb}
end
