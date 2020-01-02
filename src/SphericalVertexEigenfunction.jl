function SphericalVertexEigenfunction(domain::SphericalTriangle,
                                      vertex::Int)
    vertex >= 1 && vertex <= 3 || throw(ErrorException("vertex must be between 1 and 3 not $vertex"))
    SphericalVertexEigenfunction(domain, vertex, arb[])
end

"""
    set_eigenfunction!(u, coefficients)
> Set the coefficients for the expansion of the eigenfunction.
"""
function set_eigenfunction!(u::SphericalVertexEigenfunction,
                            coefficients::Vector)
    copy!(u.coefficients, u.domain.parent.(coefficients))
end

"""
    mu(eigenfunction::SphericalVertexEigenfunction,
       k::Integer = 1)
> Return k*μ0 as an arb ball, which is the parameter used for the
  Legendre function.
"""
function mu(u::SphericalVertexEigenfunction{fmpq},
            k::Integer = 1)
    u.domain.parent(-k*inv(u.domain.angles[u.vertex]))
end

function mu(u::SphericalVertexEigenfunction{arb},
            k::Integer = 1)
    u.domain.parent(-k*inv(u.domain.angles[u.vertex]))*u.domain.parent(π)
end

"""
    u(θ::arb, ϕ::arb, λ::arb, k::Integer)
> Evaluate the k-th basis function for the eigenfunction with the
  given λ on the point given by (θ, ϕ).
"""
function (u::SphericalVertexEigenfunction)(θ::arb,
                                           ϕ::arb,
                                           λ::arb,
                                           k::Integer)
    ν = θ.parent(-0.5) + sqrt(θ.parent(0.25) + λ)
    μ = θ.parent(mu(u, k))
    legendre_p_safe(ν, μ, cos(θ))*sin(μ*ϕ)
end

function (u::SphericalVertexEigenfunction)((θ, ϕ),
                                           λ::arb,
                                           k::Integer)
    u(θ, ϕ, λ, k)
end

"""
    u(θ::arb, ϕ::arb, λ::arb, k::Integer)
> Evaluate the eigenfunction with the given λ on the point given by
  (θ, ϕ).
"""
function (u::SphericalVertexEigenfunction)(θ::arb,
                                           ϕ::arb,
                                           λ::arb)
    res = θ.parent(0)

    for k in 1:length(u.coefficients)
        res += u.coefficients[k]*u(θ, ϕ, λ, k)
    end

    res
end

function (u::SphericalVertexEigenfunction)((θ, ϕ),
                                           λ::arb)
    u(θ, ϕ, λ)
end
