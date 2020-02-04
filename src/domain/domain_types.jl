abstract type AbstractDomain end

abstract type AbstractSphericalDomain <: AbstractDomain end

abstract type AbstractPlanarDomain <: AbstractDomain end

"""
    SphericalTriangle((α::arb, β::arb, γ::arb)[, RR::ArbField = RealField(64)])
    SphericalTriangle((α::fmpq, β::fmpq, γ::fmpq)[, RR::ArbField = RealField(64)])
> Create a spherical triangle with angles α, β, γ based on the field
  RR.

  If the angles are of type arb then α, β and γ represents the angles
  directly. If they are of type fmpq then they represent the angle as
  a rational multiple of π. The vertex with angle α is places on the
  north pole and the vertex with angle β is taken to have y = 0.
  Vertex 1, 2 and 3 correspond to angles α, β and γ respectively. The
  boundaries are ordered by which vertex they are opposite of, so
  boundary 1 is opposite of vertex 1.
"""
struct SphericalTriangle{T <: Union{fmpq, arb}} <: AbstractSphericalDomain
    angles::Tuple{T, T, T}
    parent::ArbField
end

"""
    LShape()
> Create the classical L-shaped domain.

  The only relevant vertex is the re-entrant one with angle 3π/2 and
  the boundaries are enumerated from 1 to 4 in positive direction
  skipping the two which are directly adjacent to the re-entrant
  vertex.
"""
struct LShape <: AbstractPlanarDomain
    parent::ArbField
end
