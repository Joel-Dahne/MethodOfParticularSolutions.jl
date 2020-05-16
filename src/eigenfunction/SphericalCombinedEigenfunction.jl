function SphericalCombinedEigenfunction(domain::SphericalTriangle,
                                        us::Vector{<:AbstractSphericalEigenfunction})
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

function set_domain!(u::SphericalCombinedEigenfunction,
                     domain::AbstractDomain)
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

function active_boundaries(u::SphericalCombinedEigenfunction)
    reduce((x, y) -> x .| y, active_boundaries.(u.us))
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
function active_eigenfunctions(u::SphericalCombinedEigenfunction,
                               i::Integer)
    indices = [j for j in 1:length(u.us) if active_boundaries(u.us[j])[i]]

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
function basis_function(u::SphericalCombinedEigenfunction,
                        k::Integer)
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

    i, u.orders[i]*j + k
end

function coefficients(u::SphericalCombinedEigenfunction)
    coeffs = [coefficients(v) for v in u.us]
    N = sum(length, coeffs)
    [coeffs[i][j] for (i, j) in [basis_function(u, k) for k in 1:N]]
end

function set_eigenfunction!(u::SphericalCombinedEigenfunction,
                            coefficients::Vector)
    is = [basis_function(u, k)[1] for k in 1:length(coefficients)]
    for i in 1:length(u.us)
        set_eigenfunction!(u.us[i], coefficients[is .== i])
    end
end

function (u::SphericalCombinedEigenfunction)(xyz::AbstractVector{T},
                                             λ::arb,
                                             k::Integer;
                                             boundary = nothing,
                                             notransform::Bool = false
                                             )  where {T <: Union{arb, arb_series}}
    i, j = basis_function(u, k)
    u.us[i](xyz, λ, j, boundary = boundary, notransform = notransform)
end

function (u::SphericalCombinedEigenfunction)(θ::T,
                                             ϕ::T,
                                             λ::arb,
                                             k::Integer;
                                             boundary = nothing,
                                             notransform::Bool = false
                                             ) where {T <: Union{arb, arb_series}}
    i, j = basis_function(u, k)
    u.us[i](θ, ϕ, λ, j, boundary = boundary, notransform = notransform)
end

function (u::SphericalCombinedEigenfunction)(xyz::AbstractVector{T},
                                             λ::arb;
                                             boundary = nothing,
                                             notransform::Bool = false
                                             ) where {T <: Union{arb, arb_series}}
    res = u.domain.parent(0)

    for v in u.us
        res += v(xyz, λ, boundary = boundary, notransform = notransform)
        if (T == arb && !isfinite(res)) || (T == arb_series && !isfinite(res[end]))
            return res
        end
    end

    res
end

function (u::SphericalCombinedEigenfunction)(θ::T,
                                             ϕ::T,
                                             λ::arb;
                                             boundary = nothing,
                                             notransform::Bool = false
                                             ) where {T <: Union{arb, arb_series}}
    res = u.domain.parent(0)

    for v in u.us
        res += v(θ, ϕ, λ, boundary = boundary, notransform = notransform)
        if (T == arb && !isfinite(res)) || (T == arb_series && !isfinite(res[end]))
            return res
        end
    end

    res
end
