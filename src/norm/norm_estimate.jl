function norm_estimate(
    domain::AbstractDomain,
    u::AbstractEigenfunction{S},
    λ;
    n = 1000,
    return_values = false,
    threaded = true,
) where {S}
    pts = interior_points(domain, n)
    values = similar(pts, S)

    if threaded
        Threads.@threads for i in eachindex(pts)
            values[i] = abs(u(pts[i], λ))^2
        end
    else
        for i in eachindex(pts)
            values[i] = abs(u(pts[i], λ))^2
        end
    end

    res = sqrt(area(domain) * sum(values) / n)

    if return_values
        return res, values
    else
        return res
    end
end
