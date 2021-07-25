abstract type AbstractPlanarEigenfunction <: AbstractEigenfunction end

abstract type AbstractStandalonePlanarEigenfunction{
    S<:Union{AbstractFloat,arb},
    T<:Union{AbstractFloat,arb,Rational,fmpq},
} <: AbstractPlanarEigenfunction end

struct StandaloneVertexEigenfunction{S,T} <: AbstractStandalonePlanarEigenfunction{S,T}
    vertex::SVector{2,S}
    orientation::T
    θ::T
    stride::Int
    offset::Int
    reversed::Bool
    coefficients::Vector{S}
    parent::Union{ArbField,Nothing}

    function StandaloneVertexEigenfunction(
        vertex::AbstractVector{S},
        orientation::T,
        θ::T;
        stride::Integer = 1,
        offset::Integer = 0,
        reversed::Bool = false,
        parent::Union{ArbField,Nothing} = eltype(vertex) == arb ? parent(vertex[1]) :
                                          nothing,
    ) where {S<:Union{AbstractFloat,arb},T<:Union{AbstractFloat,Rational,arb,fmpq}}
        # Prefer fmpq over Rational if S is arb
        if S == arb && T <: Rational
            U = fmpq
            orientation = fmpq(orientation)
            θ = fmpq(θ)
        else
            U = T
        end
        return new{S,U}(vertex, orientation, θ, stride, offset, reversed, [], parent)
    end

    function StandaloneVertexEigenfunction{S,T}(
        vertex::AbstractVector,
        orientation::T,
        θ::T;
        stride::Integer = 1,
        offset::Integer = 0,
        reversed::Bool = false,
        parent::Union{ArbField,Nothing} = eltype(vertex) == arb ? parent(vertex[1]) :
                                          nothing,
    ) where {S,T}
        return new{S,T}(vertex, orientation, θ, stride, offset, reversed, [], parent)
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
struct StandaloneInteriorEigenfunction{S,T} <: AbstractStandalonePlanarEigenfunction{S,T}
    vertex::SVector{2,S}
    orientation::T
    stride::Int
    offset::Int
    even::Bool
    odd::Bool
    coefficients::Vector{S}
    parent::Union{ArbField,Nothing}

    function StandaloneInteriorEigenfunction(
        vertex::AbstractVector{S},
        orientation::T = zero(ifelse(eltype(vertex) == arb, fmpq, Rational{Int}));
        offset::Integer = 0,
        stride::Integer = 1,
        even::Bool = false,
        odd::Bool = false,
        parent::Union{ArbField,Nothing} = eltype(vertex) == arb ? parent(vertex[1]) :
                                          nothing,
    ) where {S<:Union{AbstractFloat,arb},T<:Union{AbstractFloat,Rational,arb,fmpq}}
        even && odd && throw(ArgumentError("eigenfunction can't be both even and odd"))

        # Prefer fmpq over Rational if S is arb
        if S == arb && T <: Rational
            U = fmpq
            orientation = fmpq(orientation)
        else
            U = T
        end

        return new{S,U}(vertex, orientation, stride, offset, even, odd, arb[], parent)
    end

    function StandaloneInteriorEigenfunction{S,T}(
        vertex::AbstractVector,
        orientation::T = zero(T);
        offset::Integer = 0,
        stride::Integer = 1,
        even::Bool = false,
        odd::Bool = false,
        parent::Union{ArbField,Nothing} = eltype(vertex) == arb ? parent(vertex[1]) :
                                          nothing,
    ) where {S,T}
        even && odd && throw(ArgumentError("eigenfunction can't be both even and odd"))

        return new{S,T}(vertex, orientation, stride, offset, even, odd, [], parent)
    end
end

struct StandaloneLightningEigenfunction{S,T} <: AbstractStandalonePlanarEigenfunction{S,T}
    vertex::SVector{2,S}
    orientation::T
    θ::T
    l::S
    σ::S
    even::Bool
    odd::Bool
    reversed::Bool
    coefficients::Vector{S}
    parent::Union{ArbField,Nothing}

    function StandaloneLightningEigenfunction(
        vertex::AbstractVector{S},
        orientation::T,
        θ::T;
        l = 1,
        σ = 4,
        even::Bool = false,
        odd::Bool = false,
        reversed::Bool = false,
        parent::Union{ArbField,Nothing} = eltype(vertex) == arb ? parent(vertex[1]) :
                                          nothing,
    ) where {S<:Union{AbstractFloat,arb},T<:Union{AbstractFloat,Rational,arb,fmpq}}
        even && odd && throw(ArgumentError("eigenfunction can't be both even and odd"))

        # Prefer fmpq over Rational if S is arb
        if S == arb && T <: Rational
            U = fmpq
            orientation = fmpq(orientation)
            θ = fmpq(θ)
        else
            U = T
        end

        # convert(arb, x) doesn't work
        if S == arb
            l = parent(l)
            σ = parent(σ)
        end

        return new{S,U}(vertex, orientation, θ, l, σ, even, odd, reversed, [], parent)
    end

    function StandaloneLightningEigenfunction{S,T}(
        vertex::AbstractVector,
        orientation::T,
        θ::T;
        l = 1,
        σ = 4,
        even::Bool = false,
        odd::Bool = false,
        reversed::Bool = false,
        parent::Union{ArbField,Nothing} = eltype(vertex) == arb ? parent(vertex[1]) :
                                          nothing,
    ) where {S,T}
        even && odd && throw(ArgumentError("eigenfunction can't be both even and odd"))

        # convert(arb, x) doesn't work
        if S == arb
            l = parent(l)
            σ = parent(σ)
        end

        return new{S,T}(vertex, orientation, θ, l, σ, even, odd, reversed, [], parent)
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
