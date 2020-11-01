function maximize(domain::AbstractDomain,
                  u::AbstractEigenfunction,
                  λ::arb,
                  n::Integer;
                  store_trace = false,
                  show_progress = false)
    N = length(coefficients(u))
    f = t -> u(boundary_parameterization(t, domain, n), λ, boundary = n)

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
                   zero(λ),
                   one(λ),
                   absmax = true,
                   evaltype = :taylor,
                   n = N,
                   atol = 0,
                   rtol = 1e-2,
                   store_trace = store_trace,
                   extended_trace = store_trace,
                   callback = enclose_progress)
end

function maximize(domain::AbstractDomain,
                  u::AbstractEigenfunction,
                  λ::arb;
                  store_trace = false,
                  show_progress = false)
    boundaries = active_boundaries(u)
    m = zero(λ)

    traces = Dict()
    for boundary in findall(boundaries)
        m2 = maximize(domain, u, λ, boundary,
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
    @timeit_debug "maximize" m = maximize(domain, u, λ,
                                          store_trace = extended_trace,
                                          show_progress = show_progress)
    if extended_trace
        m, maximize_trace = m
    end
    @timeit_debug "norm" begin
        if norm_rigorous
            n = norm(domain, u, λ)
        else
            n = norm(domain, u, λ, numpoints = 4length(coefficients(u)))
        end
    end
    ϵ = sqrt(area(domain))*m/n
    lower = λ/(1 + getinterval(ϵ)[2])
    upper = λ/(1 - getinterval(ϵ)[2])
    if lower <= upper
        enclosure = setinterval(lower, upper)
    else
        enclosure = ball(λ, domain.parent(Inf))
    end

    if extended_trace
        return enclosure, n, m, maximize_trace
    elseif store_trace
        return enclosure, n, m
    else
       return enclosure
    end
end
