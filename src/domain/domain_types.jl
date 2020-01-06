abstract type AbstractDomain end

abstract type AbstractSphericalDomain <: AbstractDomain end

struct SphericalTriangle{T <: Union{fmpq, arb}} <: AbstractSphericalDomain
    angles::Tuple{T, T, T}
    parent::ArbField
end
