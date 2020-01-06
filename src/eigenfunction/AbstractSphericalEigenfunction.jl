"""
    set_eigenfunction!(u::AbstractSphericalEigenfunction, coefficients)
> Set the coefficients for the expansion of the eigenfunction.
"""
function set_eigenfunction!(u::AbstractSphericalEigenfunction,
                            coefficients::Vector)
    copy!(u.coefficients, u.domain.parent.(coefficients))
end

"""
    coordinate_transformation(u::AbstractVertexEigenfunction, θ::arb, ϕ::arb)
> Perform the coordinate transformation given by
  coordinate_transformation(u) directly on θ and ϕ. Depending on the
  eigenfunction this can be more efficient if the transformation is
  not used multiple times.
"""
function coordinate_transformation(u::AbstractSphericalEigenfunction,
                                   θ::arb,
                                   ϕ::arb)
    coordinate_transformation(u)(θ, ϕ)
end
