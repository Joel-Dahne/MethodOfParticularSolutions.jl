function maximize(u::AbstractEigenfunction,
                  λ::arb)
    @error "no rigorous implementation of maximize for $(typeof(u)), computing approximate maximum"
    boundary = boundary_points(u.domain, u, 1000)
    maximum(abs(u(b, λ)) for b in boundary)
end

function maximize(u::AbstractSphericalEigenfunction,
                  λ::arb,
                  i::Integer)
    N = length(coefficients(u))
    f = t -> u(boundary_parameterization(t, u.domain, i), λ)
    M = enclosemaximum(f,
                       u.domain.parent(0),
                       u.domain.parent(1),
                       absmax = true,
                       evaltype = :taylor,
                       n = N,
                       atol = 0,
                       rtol = 1e-2,
                       show_trace = false)
end

function maximize(u::AbstractSphericalEigenfunction,
                  λ::arb)
    active = active_boundaries(u)

    m = u.domain.parent(0)

    for vertex in 1:3
        if active[vertex]
            m = max(m, maximize(u, λ, vertex))
        end
    end

    m
end

function maximize(u::SphericalCombinedEigenfunction,
                  λ::arb)
    active = active_boundaries(u)

    m = u.domain.parent(0)

    for vertex in 1:3
        if active[vertex]
            v = active_eigenfunctions(u, vertex)
            m = max(m, maximize(v, λ, vertex))
        end
    end

    m
end

function enclose_eigenvalue(domain::AbstractDomain,
                            u::AbstractEigenfunction,
                            λ::arb)
    m = maximize(u, λ)
    n = norm(u, λ)
    ϵ = sqrt(area(domain))*m/n

    enclosure = λ/ball(domain.parent(1), ϵ)
end
