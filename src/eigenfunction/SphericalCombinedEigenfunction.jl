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

"""
    norm(u::SphericalCombinedEigenfunction,
         λ::arb,
         (a, b, c))
> Lower bound the norm on the spherical triangle given by the three
  vertices `a`, `b` and `c`.

  The bound is computed by lower bounding `u^2` on the boundary and
  using a version of the minimum principle. To be able to use the
  minimum principle the area of the triangle needs to be small enough
  compared to the eigenvalue.
"""
function norm(u::SphericalCombinedEigenfunction,
              λ::arb,
              (a, b, c))
    # Check that the area is small enough
    area = sum(anglesfromvertices(a, b, c)) - u.domain.parent(π)
    # FIXME: Use the correct maximum area here
    maximumarea = 4*u.domain.parent(π)^2

    if !(area < maximumarea)
        return u.domain.parent(0)
    end

    m = u.domain.parent(Inf)
    for (v, w) in [(a, b), (b, c), (c, a)]
        f = t -> -u(normalize(v .+ t.*(w - v)), λ)^2

        M = -enclosemaximum(f,
                            u.domain.parent(0),
                            u.domain.parent(1),
                            evaltype = :taylor,
                            n = 4,
                            atol = 1e-5,
                            rtol = 1e-5,
                            maxevals = 10000)

        m = min(m, M)
        if !isfinite(m) || !(m > 0)
            break
        end
    end

    res = m*area

    if isfinite(res) && res > 0
        return res
    else
        return u.domain.parent(0)
    end
end

function norm(u::SphericalCombinedEigenfunction,
              λ::arb,
              (a, b, c),
              recursions::Int)
    if recursions == 0
        res = norm(u, λ, (a, b, c))
        return res
    end

    x = normalize(a + b)
    y = normalize(b + c)
    z = normalize(c + a)

    sum(norm(u, λ, vs, recursions - 1) for vs in [(a, x, z), (x, b, y),
                                                  (x, y, z), (z, y, c)])
end

function norm(u::SphericalCombinedEigenfunction,
              λ::arb)
    n = 1
    a = normalize(boundary_points(u.domain, 2, n)[end] + boundary_points(u.domain, 3, n)[1])
    b = normalize(boundary_points(u.domain, 3, n)[end] + boundary_points(u.domain, 1, n)[1])
    c = normalize(boundary_points(u.domain, 1, n)[end] + boundary_points(u.domain, 2, n)[1])

    res = norm(u, λ, (a, b, c), 2)

    res
end
