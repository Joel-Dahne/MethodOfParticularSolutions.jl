function SphericalCombinedEigenfunction(
    domain::SphericalTriangle,
    us::Vector{<:AbstractSphericalEigenfunction},
)
    orders = ones(Int, length(us))
    !isempty(us) || throw(ArgumentError("us must be non-empty"))

    SphericalCombinedEigenfunction(domain, us, orders)
end

function Base.show(io::IO, u::SphericalCombinedEigenfunction)
    println(io, "Combined eigenfunction")
    println(io, "domain: $(u.domain)")
    println(io, "eigenfunctions:")
    recur_io = IOContext(io, :compact => true)
    for v in u.us
        print(io, "- ")
        show(recur_io, v)
    end
    print(io, "number of set coefficients: $(length(coefficients(u)))")
end

function set_domain!(u::SphericalCombinedEigenfunction, domain::AbstractDomain)
    u.domain = domain
    for v in u.us
        set_domain!(v, domain)
    end
    return u
end

function recompute!(u::SphericalCombinedEigenfunction)
    for v in u.us
        recompute!(v)
    end
    u
end

function active_boundaries(domain::SphericalTriangle, u::SphericalCombinedEigenfunction)
    if domain === u.domain
        return union([active_boundaries(domain, v) for v in u.us]...)
    else
        return 1:3
    end
end

"""
    active_eigenfunctions(u::SphericalCombinedEigenfunction,
                          i::Integer)
> Return a SphericalCombinedEigenfunction consisting of all
  eigenfunction in u which are active on boundary i.

  If there are no such eigenfunctions then nothing is returned.

  The returned eigenfunction shares coefficients with the old one, any
  changes to the first also applies to the latter and vice versa.

  See also: [`active_boundaries`](@ref)
"""
function active_eigenfunctions(u::SphericalCombinedEigenfunction, i::Integer)
    indices = [j for j in 1:length(u.us) if i ∈ active_boundaries(u.us[j])]

    if isempty(indices)
        return nothing
    end

    vs = u.us[indices]
    orders = u.orders[indices]

    SphericalCombinedEigenfunction(u.domain, vs, orders)
end

"""
    basis_function(u::SphericalCombinedEigenfunction,
                   k::Integer)
> Return the index for the eigenfunctions that u is a combination of
  that should be used and the index to use for that eigenfunction.

  This depends on the relative order between the eigenfunctions. It
  returns i, j which corresponds to calling the j-th basis-function of
  u.us[i].
"""
function basis_function(u::SphericalCombinedEigenfunction, k::Integer)
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

    i, u.orders[i] * j + k
end

"""
    basis_function(u::SphericalCombinedEigenfunction, ks::UnitRange{Int})

Return a vector with elements of type `UnitRange{Int}` where element
`i` is the range of indices that should be used for the respective
eigenfunctions, i.e. the eigenfunction `u.us[i]` should take the range
given by the `i`th element in this vector.
"""
function basis_function(u::SphericalCombinedEigenfunction, ks::UnitRange{Int})
    ks.start == 1 || throw(ArgumentError("ks must start with 1, got ks = $ks"))
    ks.stop == 0 && return fill(1:0, length(u.us))

    A = [0; cumsum(u.orders)]
    fullcycles = div(ks.stop, A[end])
    R = ks.stop % A[end]
    remaining = [max(min(u.orders[i], R - A[i]), 0) for i in eachindex(u.us)]
    return [1:(fullcycles*u.orders[i]+remaining[i]) for i in eachindex(u.us)]
end

function coefficients(u::SphericalCombinedEigenfunction)
    coeffs = [coefficients(v) for v in u.us]
    N = sum(length, coeffs)
    [coeffs[i][j] for (i, j) in [basis_function(u, k) for k = 1:N]]
end

function set_eigenfunction!(u::SphericalCombinedEigenfunction, coefficients::Vector)
    is = [basis_function(u, k)[1] for k = 1:length(coefficients)]
    for i = 1:length(u.us)
        set_eigenfunction!(u.us[i], coefficients[is.==i])
    end
end

function (u::SphericalCombinedEigenfunction)(
    xyz::AbstractVector{T},
    λ::arb,
    k::Integer;
    boundary = nothing,
    notransform::Bool = false,
) where {T<:Union{arb,arb_series}}
    i, j = basis_function(u, k)
    u.us[i](xyz, λ, j, boundary = boundary, notransform = notransform)
end

function (u::SphericalCombinedEigenfunction)(
    θ::T,
    ϕ::T,
    λ::arb,
    k::Integer;
    boundary = nothing,
    notransform::Bool = false,
) where {T<:Union{arb,arb_series}}
    i, j = basis_function(u, k)
    u.us[i](θ, ϕ, λ, j, boundary = boundary, notransform = notransform)
end

function (u::SphericalCombinedEigenfunction)(
    xyz::AbstractVector{T},
    λ::arb,
    ks::UnitRange{Int};
    boundary = nothing,
    notransform::Bool = false,
) where {T<:Union{arb,arb_series}}
    ks_per_index = basis_function(u, ks)

    res_per_index = similar(u.us, Vector{T})
    for i in eachindex(u.us)
        res_per_index[i] =
            u.us[i](xyz, λ, ks_per_index[i], boundary = boundary, notransform = notransform)
    end

    res = similar(ks, T)
    # TODO: This can be done more efficiently
    for i in eachindex(ks)
        j, l = basis_function(u, ks[i])
        res[i] = res_per_index[j][l]
    end

    return res
end

function (u::SphericalCombinedEigenfunction)(
    θ::T,
    ϕ::T,
    λ::arb,
    ks::UnitRange{Int};
    boundary = nothing,
    notransform::Bool = false,
) where {T<:Union{arb,arb_series}}
    ks_per_index = basis_function(u, ks)

    res_per_index = similar(u.us, Vector{T})
    for i in eachindex(u.us)
        res_per_index[i] = u.us[i](
            θ,
            ϕ,
            λ,
            ks_per_index[i],
            boundary = boundary,
            notransform = notransform,
        )
    end

    res = similar(ks, T)
    # TODO: This can be done more efficiently
    for i in eachindex(ks)
        j, l = basis_function(u, ks[i])
        res[i] = res_per_index[j][l]
    end

    return res
end

function (u::SphericalCombinedEigenfunction)(
    xyz::AbstractVector{T},
    λ::arb;
    boundary = nothing,
    notransform::Bool = false,
) where {T<:Union{arb,arb_series}}
    res = u.domain.parent(0)

    for v in u.us
        res += v(xyz, λ, boundary = boundary, notransform = notransform)
        if (T == arb && !isfinite(res)) || (T == arb_series && !isfinite(res[end]))
            return res
        end
    end

    res
end

function (u::SphericalCombinedEigenfunction)(
    θ::T,
    ϕ::T,
    λ::arb;
    boundary = nothing,
    notransform::Bool = false,
) where {T<:Union{arb,arb_series}}
    res = u.domain.parent(0)

    for v in u.us
        res += v(θ, ϕ, λ, boundary = boundary, notransform = notransform)
        if (T == arb && !isfinite(res)) || (T == arb_series && !isfinite(res[end]))
            return res
        end
    end

    res
end
