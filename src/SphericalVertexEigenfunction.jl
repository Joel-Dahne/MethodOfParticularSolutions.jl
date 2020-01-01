function SphericalVertexEigenfunction(domain::SphericalTriangle,
                                      vertex::Int)
    vertex >= 1 && vertex <= 3 || throw(ErrorException("vertex must be between 1 and 3 not $vertex"))
    SphericalVertexEigenfunction(domain, vertex, arb[])
end

function set_eigenfunction!(u::SphericalVertexEigenfunction,
                            coefficients::Vector)
    u.coefficients = u.domain.parent.(coefficients)
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

function (u::SphericalVertexEigenfunction)(θ::arb,
                                           ϕ::arb,
                                           λ,
                                           k::Integer)
    ν = θ.parent(-0.5) + sqrt(θ.parent(0.25) + λ)
    μ = θ.parent(mu(u, k))
    legendre_p_safe(ν, μ, cos(θ))*sin(μ*ϕ)
end

function (u::SphericalVertexEigenfunction)((θ, ϕ),
                                           λ,
                                           k::Integer)
    ν = θ.parent(-0.5) + sqrt(θ.parent(0.25) + λ)
    μ = θ.parent(mu(u, k))
    legendre_p_safe(ν, μ, cos(θ))*sin(μ*ϕ)
end

function (u::SphericalVertexEigenfunction)(θ::arb,
                                           ϕ::arb,
                                           λ)
    sum(u.coefficients .* u.(θ, ϕ, λ, 1:length(u.coefficients)))
end

function (u::SphericalVertexEigenfunction)((θ, ϕ),
                                           λ)
    res = θ.parent(0)

    for k in 1:length(u.coefficients)
        res += u.coefficients[k]*u((θ, ϕ), λ, k)
    end

    res
end
