function defect_estimate(
    domain::AbstractDomain,
    u::AbstractEigenfunction{S},
    λ;
    n = 1000,
    return_values = false,
    threaded = true,
) where {S}
    pts, idxs = boundary_points(domain, u, length(coefficients(u)), n)
    values = similar(pts, S)

    if threaded
        Threads.@threads for i in eachindex(pts)
            values[i] = u(pts[i], λ, boundary = idxs[i])
        end
    else
        for i in eachindex(pts)
            values[i] = u(pts[i], λ, boundary = idxs[i])
        end
    end

    res = maximum(abs.(values))

    if return_values
        return res, values
    else
        return res
    end
end
