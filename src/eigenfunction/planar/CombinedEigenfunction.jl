CombinedEigenfunction(
    domain::AbstractPlanarDomain{S,T},
    us::Vector{<:AbstractPlanarEigenfunction},
    orders::Vector{<:Integer} = ones(Int, length(us));
    even_boundaries = BitSet(),
    us_to_boundary = [active_boundaries(domain, u) for u in us],
) where {S,T} =
    CombinedEigenfunction{S,T}(domain, us, orders; even_boundaries, us_to_boundary)

function CombinedEigenfunction{S,T}(
    domain::AbstractPlanarDomain,
    us::Vector{<:AbstractPlanarEigenfunction},
    orders::Vector{<:Integer} = ones(Int, length(us));
    even_boundaries = BitSet(),
    us_to_boundary = [active_boundaries(domain, u) for u in us],
) where {S,T}
    boundary_to_us = OrderedDict(
        i => BitSet(findall(B -> i ∈ B, us_to_boundary)) for i in boundaries(domain)
    )

    return CombinedEigenfunction{S,T}(domain, us, orders, even_boundaries, boundary_to_us)
end


function Base.getproperty(u::CombinedEigenfunction, name::Symbol)
    if name == :parent
        return u.domain.parent
    else
        return getfield(u, name)
    end
end

function Base.show(io::IO, ::MIME"text/plain", u::CombinedEigenfunction)
    print(
        io,
        "Combined eigenfunction on $(typeof(u.domain)) with $(length(coefficients(u))) coefficients",
    )
    if !get(io, :compact, false)
        println(io, "")
        recur_io = IOContext(io, :compact => true)
        for v in u.us
            print(io, "- ")
            show(recur_io, v)
        end
    end
end

Base.show(io::IO, u::CombinedEigenfunction) = print(
    io,
    "Combined eigenfunction on $(typeof(u.domain)) with $(length(coefficients(u))) coefficients",
)

function set_domain!(u::CombinedEigenfunction, domain::AbstractDomain)
    u.domain = domain
    for v in u.us
        set_domain!(v, domain)
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

    return i, u.orders[i] * j + k
end

"""
    basis_function(u::CombinedEigenfunction, ks::UnitRange{Int})

Return a vector with elements of type `UnitRange{Int}` where element
`i` is the range of indices that should be used for the respective
eigenfunctions, i.e. the eigenfunction `u.us[i]` should take the range
given by the `i`th element in this vector.
"""
function basis_function(u::CombinedEigenfunction, ks::UnitRange{Int})
    if ks.start > 1
        res1 = basis_function(u, 1:ks.start-1)
        res2 = basis_function(u, 1:ks.stop)

        return [r1.stop+1:r2.stop for (r1, r2) in zip(res1, res2)]
    end

    ks.start == 1 ||
        throw(ArgumentError("ks must start with a positive number, got ks = $ks"))
    ks.stop == 0 && return fill(1:0, length(u.us))

    A = [0; cumsum(u.orders)]
    fullcycles = div(ks.stop, A[end])
    R = ks.stop % A[end]
    remaining = [max(min(u.orders[i], R - A[i]), 0) for i in eachindex(u.us)]
    return [1:(fullcycles*u.orders[i]+remaining[i]) for i in eachindex(u.us)]
end

function coefficients(u::CombinedEigenfunction)
    coeffs = [coefficients(v) for v in u.us]
    N = sum(length, coeffs)
    return [coeffs[i][j] for (i, j) in [basis_function(u, k) for k = 1:N]]
end

function set_eigenfunction!(u::CombinedEigenfunction, coefficients::Vector)
    is = [basis_function(u, k)[1] for k = 1:length(coefficients)]
    for i = 1:length(u.us)
        set_eigenfunction!(u.us[i], coefficients[is.==i])
    end
    return u
end

function (u::CombinedEigenfunction{S,T})(
    xy::AbstractVector,
    λ,
    ks::UnitRange{Int};
    boundary = nothing,
    notransform::Bool = false,
) where {S,T}
    # The individual eigenfunctions do not necessarily have a common
    # type and the conversion here might be redone for each
    # eigenfunction. Still we opt for an explicit conversion here to
    # give a more uniform behaviour.
    if S == arb
        if eltype(xy) == arb_series
            xy = convert(SVector{2,arb_series}, xy)
        else
            xy = convert(SVector{2,arb}, u.parent.(xy))
        end
        λ = u.parent(λ)
    else
        xy = convert(SVector{2,S}, xy)
        λ = convert(S, λ)
    end

    indices = isnothing(boundary) ? eachindex(u.us) : u.boundary_to_us[boundary]

    ks_per_index = basis_function(u, ks)
    res_per_index = similar(u.us, Vector{eltype(xy)})

    for i in eachindex(u.us)
        if i ∈ indices && !isempty(ks_per_index[i])
            res_tmp = u.us[i](xy, λ, ks_per_index[i]; boundary, notransform)
            # We don't have to handle arb_series specially here. The
            # only way res_tmp could contain arb_series is if xy
            # contains arb_series, in which case the else statement
            # works.
            if eltype(eltype(res_per_index)) == arb
                res_per_index[i] = u.parent.(res_tmp)
            else
                res_per_index[i] = res_tmp
            end
        else
            res_per_index[i] = [zero(first(xy)) for _ in eachindex(ks)]
        end
    end

    # TODO: This can be done more efficiently
    res = similar(ks, eltype(xy))
    for i in eachindex(ks)
        j, l = basis_function(u, ks[i])
        res[i] = res_per_index[j][l]
    end

    return res
end
