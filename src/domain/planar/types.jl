abstract type AbstractPlanarDomain{S,T} <: AbstractDomain{S,T} end

"""
    LShape()
> Create the classical L-shaped domain.

  The only relevant vertex is the re-entrant one with angle 3π/2 and
  the boundaries are enumerated from 1 to 4 in positive direction
  skipping the two which are directly adjacent to the re-entrant
  vertex.
"""
struct LShape <: AbstractPlanarDomain{arb,fmpq}
    parent::ArbField
end

"""
    Triangle(α, β[, parent::ArbField = RealField(64)])

Create a (planar) triangle with two of the angles being `α` and `β`.
The length of the side in between is defined to be one and the last
angle is called `γ`.

If the angles are of type `arb` then `α` and `β` represents the angles
directly. If they are of type `fmpq` then they represent the angle as
a rational multiple of `π`. The vertex with angle `α` is placed at the
origin the vertex with angle `β` is taken to have `y = 0`.

Vertex 1, 2 and 3 are opposite of angles α, β and γ respectively.
"""
struct Triangle{S,T} <: AbstractPlanarDomain{S,T}
    angles::NTuple{3,T}
    parent::Union{ArbField,Nothing}

    function Triangle(α::T, β::T; parent::Nothing = nothing) where {T<:AbstractFloat}
        α > 0 || throw(DomainError(α, "angle must be positive"))
        β > 0 || throw(DomainError(β, "angle must be positive"))
        α + β < π || throw(ArgumentError("α + β must be less than π"))
        γ = π - α - β

        return new{T,T}((α, β, γ), parent)
    end

    function Triangle(α::arb, β::arb; parent::ArbField = parent(α))
        α > 0 || throw(DomainError(α, "angle must be positive"))
        β > 0 || throw(DomainError(β, "angle must be positive"))
        !(α + β > parent(π)) || throw(ArgumentError("α + β must be less than π"))
        α + β < parent(π) || @warn "α + β might not be less than π"

        γ = parent(π) - α - β

        return new{arb,arb}((α, β, γ), parent)
    end

    function Triangle(α::T, β::T; parent::Nothing = nothing) where {T<:Rational}
        α > 0 || throw(DomainError(α, "angle must be positive"))
        β > 0 || throw(DomainError(β, "angle must be positive"))
        α + β < 1 || throw(ArgumentError("α + β must be less than 1"))
        γ = one(α) - α - β

        return new{float(T),T}((α, β, γ), parent)
    end

    function Triangle(α::fmpq, β::fmpq; parent::ArbField = RealField(64))
        α > 0 || throw(DomainError(α, "angle must be positive"))
        β > 0 || throw(DomainError(β, "angle must be positive"))
        α + β < 1 || throw(ArgumentError("α + β must be less than 1"))
        γ = 1 - α - β

        return new{arb,fmpq}((α, β, γ), parent)
    end
end

"""
    Polygon(angles, vertices)

Create a (planar) polygon with the given angles. To simplify the
implementation we also require that the location. of the vertices are
given.

The boundaries are enumerated by which vertex they are next to, in
positive order.
"""
struct Polygon{S,T} <: AbstractPlanarDomain{S,T}
    angles::Vector{T}
    vertices::Vector{SVector{2,S}}
    parent::Union{ArbField,Nothing}

    function Polygon{S,T}(
        angles::AbstractVector,
        vertices::AbstractVector;
        parent::Union{ArbField,Nothing} = nothing,
    ) where {S,T}
        length(angles) == length(vertices) || throw(
            ArgumentError(
                "length of angles and vertices don't match, got $length(angles) and length(vertices)",
            ),
        )

        angles = convert(Vector{T}, angles)
        vertices = convert(Vector{SVector{2,S}}, vertices)

        if T == arb && isnothing(parent)
            parent = Nemo.parent(vertices[1][1])
        end

        return new{S,T}(angles, vertices, parent)
    end

    function Polygon(
        angles::AbstractVector,
        vertices::AbstractVector;
        parent::Union{ArbField,Nothing} = nothing,
    )
        S = eltype(eltype(vertices))
        T = eltype(angles)

        if S != arb
            S = float(S)
        end

        return Polygon{S,T}(angles, vertices; parent)
    end
end

"""
    TransformedDomain{T,U<:AbstractPlanarDomain}(domain::U, rotation, scaling, translation)

Represents a transformation of `domain` corresponding of a rotation,
a uniform scaling and a translation (in that order).
"""
struct TransformedDomain{T<:Union{fmpq,arb},U<:AbstractPlanarDomain} <:
       AbstractPlanarDomain{arb,T}
    original::U
    rotation::T
    scaling::arb
    translation::SVector{2,arb}
end

"""
    IntersectedDomain{U<:AbstractPlanarDomain}(exterior::T, interiors::Vector{AbstractPlanarDomain})

Represents a domain given by removing from `exterior` the domains in
`interiors`. It's assumed that all domains in `interiors` are
contained in `exterior` and that they are pairwise disjoint.

The boundaries are enumerated in a way such that first comes the
exterior domains boundaries and then comes the interior domains.
"""
struct IntersectedDomain{S,T} <: AbstractPlanarDomain{S,T}
    exterior::AbstractPlanarDomain{S,T}
    interiors::Vector{<:AbstractPlanarDomain{S,T}}
end
