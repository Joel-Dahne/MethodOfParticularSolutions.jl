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
                                             notransform::Bool = false
                                             )  where {T <: Union{arb, arb_series}}
    i, j = basis_function(u, k)
    u.us[i](xyz, λ, j, notransform = notransform)
end

function (u::SphericalCombinedEigenfunction)(θ::T,
                                             ϕ::T,
                                             λ::arb,
                                             k::Integer;
                                             notransform::Bool = false
                                             ) where {T <: Union{arb, arb_series}}
    i, j = basis_function(u, k)
    u.us[i](θ, ϕ, λ, j, notransform = notransform)
end

function (u::SphericalCombinedEigenfunction)(xyz::AbstractVector{T},
                                             λ::arb;
                                             notransform::Bool = false
                                             ) where {T <: Union{arb, arb_series}}
    res = u.domain.parent(0)

    for v in u.us
        res += v(xyz, λ, notransform = notransform)
    end

    res
end

function (u::SphericalCombinedEigenfunction)(θ::T,
                                             ϕ::T,
                                             λ::arb;
                                             notransform::Bool = false
                                             ) where {T <: Union{arb, arb_series}}
    res = u.domain.parent(0)

    for v in u.us
        res += v(θ, ϕ, λ, notransform = notransform)
    end

    res
end
