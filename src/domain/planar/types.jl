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
    Polygon{S,T}(vertices, angles, orientation; parent)
    Polygon(vertices[, angles, orientation]; parent)

Create a polygon with the given vertices.

The vertices are assumed to be ordered such that when going along the
boundary from the first to the second vertex the interior of the
domain is to the left. No check is made to see if the vertices define
a proper polygon.

In addition to the vertices the polygon also stores the angles of all
the vertices as well as the orientation of the polygon (given by the
orientation of the first vertex as defined in [`orientation`](@ref)).
If the angles are stored as rational numbers these angles can
typically not be computed exactly from the vertices and they have to
be given explicitly.

The boundaries are enumerated by which two vertices they are in
between.
"""
struct Polygon{S,T} <: AbstractPlanarDomain{S,T}
    vertices::Vector{SVector{2,S}}
    angles::Vector{T}
    orientation::T
    parent::Union{ArbField,Nothing}

    function Polygon{S,T}(
        vertices::AbstractVector,
        angles::AbstractVector,
        orientation;
        parent::Union{ArbField,Nothing} = nothing,
    ) where {S,T}
        length(angles) == length(vertices) || throw(
            ArgumentError(
                "length of angles and vertices don't match, got $length(angles) and length(vertices)",
            ),
        )

        vertices = convert(Vector{SVector{2,S}}, vertices)
        angles = convert(Vector{T}, angles)
        orientation = convert(T, orientation)

        if T == arb && isnothing(parent)
            parent = Nemo.parent(vertices[1][1])
        end

        return new{S,T}(vertices, angles, orientation, parent)
    end

    function Polygon(
        vertices::AbstractVector{<:AbstractVector{S}},
        angles::AbstractVector{T},
        orientation::T;
        parent::Union{ArbField,Nothing} = nothing,
    ) where {S,T}
        length(angles) == length(vertices) || throw(
            ArgumentError(
                "length of angles and vertices don't match, got $length(angles) and length(vertices)",
            ),
        )

        vertices = convert(Vector{SVector{2,S}}, vertices)
        angles = convert(Vector{T}, angles)

        if T == arb && isnothing(parent)
            parent = Nemo.parent(vertices[1][1])
        end

        return new{S,T}(vertices, angles, orientation, parent)
    end
end

"""
    TransformedDomain{S,T,U<:AbstractPlanarDomain{S,T}}(domain::U, rotation, scaling, translation)

Represents a transformation of `domain` corresponding of a rotation,
a uniform scaling and a translation (in that order).
"""
struct TransformedDomain{
    S<:Union{AbstractFloat,arb},
    T<:Union{AbstractFloat,arb,Rational,fmpq},
    U<:AbstractPlanarDomain{S,T},
} <: AbstractPlanarDomain{S,T}
    original::U
    rotation::T
    scaling::S
    translation::SVector{2,S}

    function TransformedDomain(
        domain::AbstractPlanarDomain{S,T},
        rotation,
        scaling,
        translation,
    ) where {S,T}
        if S == arb
            scaling = domain.parent(scaling)
            translation = convert(SVector{2,S}, domain.parent.(translation))
        else
            scaling = convert(S, scaling)
            translation = convert(SVector{2,S}, translation)
        end
        if T == arb
            rotation = domain.parent(rotation)
        elseif T == fmpq
            rotation = fmpq(rotation)
        else
            # TODO: Should we do something for rational input but
            # float output and vice versa?
            rotation = convert(T, rotation)
        end

        return new{S,T,typeof(domain)}(domain, rotation, scaling, translation)
    end
end

"""
    IntersectedDomain{S,T,U<:AbstractPlanarDomain{S,T}}(exterior::U, interiors::Vector{AbstractPlanarDomain})

Represents a domain given by removing from `exterior` the domains in
`interiors`. It's assumed that all domains in `interiors` are
contained in `exterior` and that they are pairwise disjoint.

The boundaries are enumerated in a way such that first comes the
exterior domains boundaries and then comes the interior domains.
"""
struct IntersectedDomain{S,T,U<:AbstractPlanarDomain{S,T},V<:AbstractPlanarDomain{S,T}} <:
       AbstractPlanarDomain{S,T}
    exterior::U
    interiors::Vector{V}
end
