abstract type AbstractPlanarEigenfunction{S,T} <: AbstractEigenfunction{S,T} end

abstract type AbstractStandalonePlanarEigenfunction{S,T} <: AbstractPlanarEigenfunction{S,T} end

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

struct LinkedEigenfunction{S,T,U<:AbstractPlanarEigenfunction{S,T}} <:
       AbstractPlanarEigenfunction{S,T}
    us::Vector{U}
    extra_coefficients::Vector{S}
    excluded_boundaries::Vector{BitSet}

    LinkedEigenfunction(
        us::AbstractVector{U},
        extra_coefficients::AbstractVector = ones(length(us));
        excluded_boundaries::AbstractVector = fill(BitSet(), length(us)),
    ) where {S,T,U<:AbstractPlanarEigenfunction{S,T}} =
        LinkedEigenfunction{S,T,U}(us, extra_coefficients; excluded_boundaries)


    function LinkedEigenfunction{S,T,U}(
        us::AbstractVector{U},
        extra_coefficients::AbstractVector = ones(length(us));
        excluded_boundaries::AbstractVector = fill(BitSet(), length(us)),
    ) where {S,T,U}
        isempty(us) && throw(ArgumentError("us must not be empty"))
        length(us) == length(extra_coefficients) == length(excluded_boundaries) || throw(
            ArgumentError(
                "us, extra_coefficients and extra_boundaries must all have the same size",
            ),
        )

        if S == arb
            extra_coefficients = first(us).parent.(extra_coefficients)
        end

        return new{S,T,U}(us, extra_coefficients, excluded_boundaries)
    end
end

mutable struct CombinedEigenfunction{S,T} <: AbstractPlanarEigenfunction{S,T}
    domain::AbstractPlanarDomain
    us::Vector{<:AbstractPlanarEigenfunction}
    orders::Vector{Int}
    even_boundaries::BitSet
    boundary_to_us::OrderedDict{Int,BitSet}

    function CombinedEigenfunction{S,T}(
        domain,
        us,
        orders,
        even_boundaries,
        boundary_to_us,
    ) where {S,T}
        length(us) == length(orders) ||
            throw(ArgumentError("us and orders must have the same length"))

        Set(boundaries(domain)) == Set(keys(boundary_to_us)) || throw(
            ArgumentError(
                "boundary_to_us needs to have as keys the boundaries of the domain",
            ),
        )

        return new{S,T}(domain, us, orders, BitSet(even_boundaries), boundary_to_us)
    end
end

mutable struct LShapeEigenfunction <: AbstractPlanarEigenfunction{arb,fmpq}
    domain::LShape
    stride::Int
    coefficients::Vector{arb}
end
