function maximize(u::AbstractEigenfunction,
                  λ::arb;
                  store_trace = false,
                  show_trace = false,
                  extended_trace = false)
    @error "no rigorous implementation of maximize for $(typeof(u)), computing approximate maximum"
    boundary = boundary_points(u.domain, u, 1000)
    maximum(abs(u(b, λ)) for b in boundary)
end

function maximize(u::AbstractSphericalEigenfunction,
                  λ::arb,
                  i::Integer;
                  store_trace = false,
                  show_trace = false,
                  extended_trace = false)
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
                       store_trace = store_trace,
                       show_trace = show_trace,
                       extended_trace = extended_trace)
end

function maximize(u::AbstractSphericalEigenfunction,
                  λ::arb;
                  kwargs...)
    active = active_boundaries(u)

    m = u.domain.parent(0)

    for vertex in 1:3
        if active[vertex]
            m = max(m, maximize(u, λ, vertex; kwargs...))
        end
    end

    m
end

function maximize(u::SphericalCombinedEigenfunction,
                  λ::arb;
                  kwargs...)
    active = active_boundaries(u)

    m = u.domain.parent(0)

    for vertex in 1:3
        if active[vertex]
            v = active_eigenfunctions(u, vertex)
            m = max(m, maximize(v, λ, vertex; kwargs...))
        end
    end

    m
end

function enclose_eigenvalue(domain::AbstractDomain,
                            u::AbstractEigenfunction,
                            λ::arb;
                            kwargs...)
    @timeit_debug "maximize" m = maximize(u, λ; kwargs...)
    @timeit_debug "norm" n = norm(u, λ)
    ϵ = sqrt(area(domain))*m/n

    enclosure = λ/ball(domain.parent(1), ϵ)
end
