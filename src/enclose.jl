function maximize(
    domain::AbstractDomain,
    u::AbstractEigenfunction,
    λ::arb,
    n::Integer;
    store_trace = false,
    show_progress = false,
)
    N = length(coefficients(u))
    f = t -> u(boundary_parameterization(t, domain, n), λ, boundary = n)

    if show_progress
        enclose_progress(state) = begin
            sticky =
                radius(state.maxenclosure) / abs(state.maxenclosure) <= 1e-2 ? :done : true
            @info "Computing enclosure of maximum" state sticky = sticky
            false
        end
    else
        enclose_progress = nothing
    end

    enclosemaximum(
        f,
        zero(λ),
        one(λ),
        absmax = true,
        evaltype = :taylor,
        n = N,
        atol = 0,
        rtol = 1e-2,
        store_trace = store_trace,
        extended_trace = store_trace,
        callback = enclose_progress,
    )
end

function maximize(
    domain::AbstractDomain,
    u::AbstractEigenfunction,
    λ::arb;
    store_trace = false,
    show_progress = false,
)
    boundaries = active_boundaries(domain, u)
    m = zero(λ)

    traces = Dict()
    for boundary in boundaries
        m2 = maximize(
            domain,
            u,
            λ,
            boundary,
            store_trace = store_trace,
            show_progress = show_progress,
        )
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

function enclose_eigenvalue(
    domain::AbstractDomain,
    u::AbstractEigenfunction,
    λ::arb;
    rigorous_norm = true,
    show_progress = false,
    store_trace = false,
    extended_trace = false,
)
    @timeit_debug "maximize" m =
        maximize(domain, u, λ, store_trace = extended_trace, show_progress = show_progress)
    if extended_trace
        m, trace = m
    end

    @timeit_debug "norm" begin
        if rigorous_norm
            n = norm(domain, u, λ)
        else
            n = norm(domain, u, λ, numpoints = 4length(coefficients(u)), warn = false)
        end
    end

    ϵ = sqrt(area(domain)) * m / n
    lower = λ / (1 + getinterval(ϵ)[2])
    upper = λ / (1 - getinterval(ϵ)[2])
    if lower <= upper
        enclosure = setinterval(lower, upper)
    else
        enclosure = ball(λ, domain.parent(Inf))
    end

    if extended_trace
        return enclosure, n, m, trace
    elseif store_trace
        return enclosure, n, m
    else
        return enclosure
    end
end

function enclose_eigenvalue_approx(
    domain::AbstractDomain,
    u::AbstractEigenfunction,
    λ::arb;
    max_numpoints::Integer = 16length(coefficients(u)),
    norm_numpoints::Integer = 8length(coefficients(u)),
    store_trace::Bool = false,
    extended_trace::Bool = false,
)
    ## Approximate maximum on boundary
    # Points to evaluate u on
    pts, bds = boundary_points(domain, u, length(coefficients(u)), max_numpoints)
    values = similar(pts, arb)
    Threads.@threads for i in eachindex(pts)
        values[i] = u(pts[i], λ, boundary = bds[i])
    end
    m = zero(λ)
    for v in values
        m = max(m, abs(v))
    end

    ## Approximate norm
    n = norm(domain, u, λ, numpoints = norm_numpoints, warn = false)

    ## Compute enclosure
    ϵ = sqrt(area(domain)) * m / n
    lower = λ / (1 + getinterval(ϵ)[2])
    upper = λ / (1 - getinterval(ϵ)[2])
    if lower <= upper
        enclosure = setinterval(lower, upper)
    else
        enclosure = ball(λ, domain.parent(Inf))
    end

    if extended_trace
        return enclosure, n, m, Dict("values" => values, "boundaries" => bds)
    elseif store_trace
        return enclosure, n, m
    else
        return enclosure
    end
end
