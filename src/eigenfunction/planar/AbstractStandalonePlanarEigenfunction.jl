function set_eigenfunction!(
    u::AbstractStandalonePlanarEigenfunction{arb},
    coefficients::Vector,
)
    resize!(u.coefficients, length(coefficients))
    copy!(u.coefficients, u.parent.(coefficients))
    return u
end

function set_eigenfunction!(u::AbstractStandalonePlanarEigenfunction, coefficients::Vector)
    resize!(u.coefficients, length(coefficients))
    copy!(u.coefficients, coefficients)
    return u
end
