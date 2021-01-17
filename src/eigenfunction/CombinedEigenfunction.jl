function CombinedEigenfunction(
    domain::AbstractPlanarDomain,
    us::Vector{<:AbstractPlanarEigenfunction},
    orders::Vector{Int} = ones(Int, length(us));
    us_to_boundary = [active_boundaries(domain, u) for u in us],
    even_boundaries = BitSet(),
)
    boundary_to_us = OrderedDict(
        i => BitSet(findall(B -> i ∈ B, us_to_boundary))
        for i in boundaries(domain)
    )

    return CombinedEigenfunction(
        domain,
        us,
        boundary_to_us,
        BitSet(even_boundaries),
        orders,
    )
end

function CombinedEigenfunction(
    domain::IntersectedDomain,
    orders::Vector{Int} = ones(Int, length(boundaries(domain))),
)
    boundary_to_us = OrderedDict(i => BitSet() for i in boundaries(domain))
    us = AbstractPlanarEigenfunction[]

    for i in exterior_boundaries(domain)
        # TODO: Handle other types
        @assert typeof(domain.exterior) <: Triangle ||
            typeof(domain.exterior) <: TransformedDomain{<:Triangle}

        v = StandaloneVertexEigenfunction(domain.exterior, i)
        push!(us, v)
        idx = length(us)
        push!(boundary_to_us[i], idx)
        for j in interior_boundaries(domain)
            push!(boundary_to_us[j], idx)
        end
    end

    for i in interior_boundaries(domain)
        _, ip = get_domain_and_boundary(domain, i)
        # TODO: Handle other types
        @assert typeof(domain.exterior) <: Triangle ||
            typeof(domain.exterior) <: TransformedDomain{<:Triangle}

        v = StandaloneVertexEigenfunction(domain.interior, ip, outside = true)
        push!(us, v)

        idx = length(us)
        push!(boundary_to_us[i], length(us))
        for j in exterior_boundaries(domain)
            push!(boundary_to_us[j], idx)
        end
    end

    # TODO: Add interior eigenfunction
    v = StandaloneInteriorEigenfunction(domain)
    push!(us, v)
    idx = length(us)
    for j in boundaries(domain)
        push!(boundary_to_us[j], idx)
    end

    return CombinedEigenfunction(
        domain,
        us,
        boundary_to_us,
        BitSet(),
        orders,
    )
end

function Base.show(io::IO, u::CombinedEigenfunction)
    println(io, "Combined eigenfunction on $(u.domain)")
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

"""
    basis_function(u::CombinedEigenfunction, ks::UnitRange{Int})

Return a vector with elements of type `UnitRange{Int}` where element
`i` is the range of indices that should be used for the respective
eigenfunctions, i.e. the eigenfunction `u.us[i]` should take the range
given by the `i`th element in this vector.
"""
function basis_function(u::CombinedEigenfunction, ks::UnitRange{Int})
    ks.start == 1 || throw(ArgumentError("ks must start with 1, got ks = $ks"))
    ks.stop == 0 && return fill(1:0, length(u.us))

    A = [0; cumsum(u.orders)]
    fullcycles = div(ks.stop, A[end])
    R = ks.stop%A[end]
    remaining = [max(min(u.orders[i], R - A[i]), 0) for i in eachindex(u.us)]
    return [
        1:(fullcycles*u.orders[i] + remaining[i])
        for i in eachindex(u.us)
    ]
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

function (u::CombinedEigenfunction)(xy::AbstractVector{T},
                                    λ::arb,
                                    ks::UnitRange{Int};
                                    boundary = nothing,
                                    notransform::Bool = false
                                    ) where {T <: Union{arb, arb_series}}
    if isnothing(boundary)
        indices = eachindex(u.us)
    else
        indices = u.boundary_to_us[boundary]
    end

    ks_per_index = basis_function(u, ks)

    res_per_index = similar(u.us, Vector{T})
    for i in eachindex(u.us)
        if i in indices
            res_per_index[i] = u.us[i](
                xy,
                λ,
                ks_per_index[i],
                boundary = boundary,
                notransform = notransform,
            )
        else
            z = T == arb ? zero(λ) : 0*xy[1]
            res_per_index[i] = fill(z, length(ks_per_index[i]))
        end
    end

    # TODO: This can be done more efficiently
    res = similar(ks, T)
    for i in eachindex(ks)
        j, l = basis_function(u, ks[i])
        res[i] = res_per_index[j][l]
    end

    return res
end

function (u::CombinedEigenfunction)(xy::AbstractVector{T},
                                    λ::arb;
                                    boundary = nothing,
                                    notransform::Bool = false
                                    ) where {T <: Union{arb, arb_series}}
    if isnothing(boundary)
        indices = eachindex(u.us)
    else
        indices = u.boundary_to_us[boundary]
    end

    res = u.domain.parent(0)
    for i in indices
        res += u.us[i](xy, λ, boundary = boundary, notransform = notransform)
        if (T == arb && !isfinite(res)) || (T == arb_series && !isfinite(res[end]))
            return res
        end
    end

    return res
end
