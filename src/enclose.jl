function maximize(u::AbstractEigenfunction,
                  λ::arb;
                  kwargs...)
    @error "no rigorous implementation of maximize for $(typeof(u)), computing approximate maximum"
    boundary = boundary_points(u.domain, u, 1000)
    maximum(abs(u(b, λ)) for b in boundary)
end

function maximize(u::AbstractSphericalEigenfunction,
                  λ::arb,
                  i::Integer;
                  store_trace = false,
                  show_trace = false,
                  show_progress = false,
                  extended_trace = false)
    N = length(coefficients(u))
    f = t -> u(boundary_parameterization(t, u.domain, i), λ)

    if show_progress
        enclose_progress(state) = begin
            sticky = radius(state.maxenclosure)/abs(state.maxenclosure) <= 1e-2 ? :done : true
            @info "Computing enclosure of maximum" state sticky = sticky
            false
        end
    else
        enclose_progress = nothing
    end

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
                       callback = enclose_progress,
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

function maximize(u::LShapeEigenfunction,
                  λ::arb;
                  store_trace = false,
                  show_trace = false,
                  show_progress = false,
                  extended_trace = false)
    N = length(coefficients(u))
    m = u.domain.parent(0)

    for i in 1:4
        f = t -> u(boundary_parameterization(t, u.domain, i), λ)

        m = max(m,
                enclosemaximum(f,
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
                )
    end

    m
end

function enclose_eigenvalue(domain::AbstractDomain,
                            u::AbstractEigenfunction,
                            λ::arb;
                            norm_rigorous = true,
                            kwargs...)
    @timeit_debug "maximize" m = maximize(u, λ; kwargs...)
    @timeit_debug "norm" begin
        if norm_rigorous
            n = norm(u, λ)
        else
            n = norm(u, λ, numpoints = 4length(coefficients(u)))
        end
    end
    ϵ = sqrt(area(domain))*m/n

    enclosure = λ/ball(domain.parent(1), ϵ)
end
