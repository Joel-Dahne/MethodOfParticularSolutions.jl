abstract type AbstractDomain end

abstract type AbstractSphericalDomain <: AbstractDomain end


"""
    SphericalTriangle((α::arb, β::arb, γ::arb)[, RR::ArbField = RealField(64)])
    SphericalTriangle((α::fmpq, β::fmpq, γ::fmpq)[, RR::ArbField = RealField(64)])
> Create a spherical triangle with angles α, β, γ based on the field
  RR.

  If the angles are of type arb then α, β and γ represents the angles
  directly. If they are of type fmpq then they represent the angle as
  a rational multiple of π. The vertex with angle α is places on the
  north pole and the vertex with angle β is taken to have y = 0.
"""
struct SphericalTriangle{T <: Union{fmpq, arb}} <: AbstractSphericalDomain
    angles::Tuple{T, T, T}
    parent::ArbField
end
