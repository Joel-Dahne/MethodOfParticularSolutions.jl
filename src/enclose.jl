function maximize(u::AbstractEigenfunction,
                  λ::arb,
                  n::Integer;
                  store_trace = false,
                  show_progress = false)
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

    enclosemaximum(f,
                   u.domain.parent(0),
                   u.domain.parent(1),
                   absmax = true,
                   evaltype = :taylor,
                   n = N,
                   atol = 0,
                   rtol = 1e-2,
                   store_trace = store_trace,
                   extended_trace = store_trace,
                   callback = enclose_progress)
end

function maximize(u::AbstractEigenfunction,
                  λ::arb;
                  store_trace = false,
                  show_progress = false)
    boundaries = active_boundaries(u)
    m = u.domain.parent(0)

    traces = Dict()
    for boundary in findall(boundaries)
        m2 = maximize(u, λ, boundary,
                      store_trace = store_trace,
                      show_progress = show_progress)
        if store_trace
            m2, trace = m2
            traces[boundary] = trace
        end
        m = max(m, m2)
        if !isfinite(m)
            return m
        end
    end

    if store_trace
        return m, traces
    else
        return m
    end
end

function enclose_eigenvalue(domain::AbstractDomain,
                            u::AbstractEigenfunction,
                            λ::arb;
                            norm_rigorous = true,
                            store_trace = false,
                            extended_trace = false,
                            show_progress = false)
    @timeit_debug "maximize" m = maximize(u, λ,
                                          store_trace = extended_trace,
                                          show_progress = show_progress)
    if extended_trace
        m, maximize_trace = m
    end
    @timeit_debug "norm" begin
        if norm_rigorous
            n = norm(u, λ)
        else
            n = norm(u, λ, numpoints = 4length(coefficients(u)))
        end
    end
    ϵ = sqrt(area(domain))*m/n

    enclosure = λ/ball(domain.parent(1), ϵ)

    if extended_trace
        return enclosure, n, m, maximize_trace
    elseif store_trace
        return enclosure, n, m
    else
       return enclosure
    end
end
