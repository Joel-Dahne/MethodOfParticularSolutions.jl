function maximize(u::AbstractEigenfunction,
                  λ::arb;
                  kwargs...)
    @error "no rigorous implementation of maximize for $(typeof(u)), computing approximate maximum"
    boundary = boundary_points(u.domain, u, 1000)
    maximum(abs(u(b, λ)) for b in boundary)
end

function maximize(u::AbstractEigenfunction,
                  λ::arb,
                  n::Integer;
                  store_trace = false,
                  show_trace = false,
                  show_progress = false,
                  extended_trace = false)
    N = length(coefficients(u))
    f = t -> u(boundary_parameterization(t, u.domain, n), λ, boundary = n)

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

function maximize(u::AbstractEigenfunction,
                  λ::arb;
                  kwargs...)
    boundaries = active_boundaries(u)
    m = u.domain.parent(0)

    for boundary in findall(boundaries)
        m = max(m, maximize(u, λ, boundary; kwargs...))
        if !isfinite(m)
            return m
        end
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
