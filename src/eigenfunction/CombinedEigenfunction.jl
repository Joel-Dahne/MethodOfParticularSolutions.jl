function CombinedEigenfunction(
    domain::AbstractPlanarDomain,
    us::Vector{<:AbstractPlanarEigenfunction},
    orders::Vector{Int} = ones(Int, length(us)),
)
    actbnd = [active_boundaries(domain, u) for u in us]
    boundary_to_us = OrderedDict(
        i => BitSet(findall(B -> i ∈ B, actbnd))
        for i in boundaries(domain)
    )
    return CombinedEigenfunction(
        domain,
        us,
        boundary_to_us,
        orders,
    )
end

function Base.show(io::IO, u::CombinedEigenfunction)
    println(io, "Combined eigenfunction on $(u.domain)")
    println(io, "eigenfunctions:")
    recur_io = IOContext(io, :compact => true)
    for v in u.us
        print(io, "- ")
        show(recur_io, v)
    end
    print(io, "number of set coefficients: $(length(coefficients(u)))")
end

function set_domain!(u::CombinedEigenfunction, domain::AbstractDomain)
    u.domain = domain
    for v in u.us
        set_domain!(v, domain)
    end
    return u
end

function recompute!(u::CombinedEigenfunction)
    for v in u.us
        recompute!(v)
    end
    return u
end

function active_boundaries(domain::AbstractPlanarDomain, u::CombinedEigenfunction)
    if domain === u.domain
        return findall(!isempty, u.boundary_to_us)
    else
        return boundaries(domain)
    end
end

"""
    basis_function(u::CombinedEigenfunction, k::Integer)

Return the index for the eigenfunctions that `u` is a combination of
that should be used and the index to use for that eigenfunction.

This depends on the relative order between the eigenfunctions. It
returns `(i, j)` which corresponds to calling the `j`-th
basis-function of `u.us[i]`.

TODO: This could be optimized a lot more.
"""
function basis_function(u::CombinedEigenfunction, k::Integer)
    k > 0 || throw(ArgumentError("k must be positive not $k"))
    i = 1
    j = 0
    while k > u.orders[i]
        k -= u.orders[i]

        i += 1
        if i > length(u.orders)
            i = 1
            j += 1
        end
    end

    return i, u.orders[i]*j + k
end

function coefficients(u::CombinedEigenfunction)
    coeffs = [coefficients(v) for v in u.us]
    N = sum(length, coeffs)
    return [coeffs[i][j] for (i, j) in [basis_function(u, k) for k in 1:N]]
end

function set_eigenfunction!(u::CombinedEigenfunction,
                            coefficients::Vector)
    is = [basis_function(u, k)[1] for k in 1:length(coefficients)]
    for i in 1:length(u.us)
        set_eigenfunction!(u.us[i], coefficients[is .== i])
    end
    return u
end

function (u::CombinedEigenfunction)(xy::AbstractVector{T},
                                    λ::arb,
                                    k::Integer;
                                    boundary = nothing,
                                    notransform::Bool = false
                                    )  where {T <: Union{arb, arb_series}}
    i, j = basis_function(u, k)
    if !isnothing(boundary) && !(i ∈ u.boundary_to_us[boundary])
        if T == arb
            return u.domain.parent(0)
        else
            return 0*xy[1]
        end
    end

    return u.us[i](xy, λ, j, boundary = boundary, notransform = notransform)
end

function (u::CombinedEigenfunction)(r::T,
                                    θ::T,
                                    λ::arb,
                                    k::Integer;
                                    boundary = nothing,
                                    notransform::Bool = false
                                    ) where {T <: Union{arb, arb_series}}
    i, j = basis_function(u, k)
    if !isnothing(boundary) && !(i ∈ u.boundary_to_us[boundary])
        if T == arb
            return u.domain.parent(0)
        else
            return 0*xy[1]
        end
    end

    return u.us[i](r, θ, λ, j, boundary = boundary, notransform = notransform)
end

function (u::CombinedEigenfunction)(xy::AbstractVector{T},
                                    λ::arb;
                                    boundary = nothing,
                                    notransform::Bool = false
                                    ) where {T <: Union{arb, arb_series}}
    res = u.domain.parent(0)

    if isnothing(boundary)
        indices = eachindex(u.us)
    else
        indices = u.boundary_to_us[boundary]
    end

    for i in indices
        res += u.us[i](xy, λ, boundary = boundary, notransform = notransform)
        if (T == arb && !isfinite(res)) || (T == arb_series && !isfinite(res[end]))
            return res
        end
    end

    return res
end

function (u::CombinedEigenfunction)(r::T,
                                    θ::T,
                                    λ::arb;
                                    boundary = nothing,
                                    notransform::Bool = false
                                    ) where {T <: Union{arb, arb_series}}
    res = u.domain.parent(0)

    if isnothing(boundary)
        indices = eachindex(u.us)
    else
        indices = u.boundary_to_us[boundary]
    end

    for i in indices
        res += u.us[i](r, θ, λ, boundary = boundary, notransform = notransform)
        if (T == arb && !isfinite(res)) || (T == arb_series && !isfinite(res[end]))
            return res
        end
    end

    return res
end