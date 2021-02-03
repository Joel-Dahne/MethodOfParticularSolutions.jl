"""
    certifyindices(domain::SphericalTriangle{T},
                   λs::AbstractVector{arb}) where {T}

Attempt to determine the indices for the eigenvalues `λs` of the domain.
Returns a list consisting of tuples `(λ, minindex, maxindex)` where
`minindex` and `maxindex` are the minimum and maximum possible indices
for the eigenvalue λ. In particular if `minindex = maxindex` then that
is the exact index of λ.

The `minindex` is determined simply by determining how many
eigenvalues in `λs` are smaller then `λ`. For the upper bound it
computes eigenvalues of enclosing domains and uses those to lower
bound the eigenvalues of the domain.
"""
function certifyindices(
    domain::SphericalTriangle{T},
    λs::AbstractVector{arb};
    show_trace = false,
) where {T}
    indices = certifyindices(domain, λs, 1, show_trace = show_trace)

    if all(((λ, minindex, maxindex),) -> minindex == maxindex, indices)
        return indices
    end

    for vertex = 2:3
        indices2 = certifyindices(domain, λs, vertex, show_trace = show_trace)

        indices = [
            (λs[i], max(indices[i][2], indices2[i][2]), min(indices[i][3], indices2[i][3])) for i = 1:length(indices)
        ]

        if all(((λ, minindex, maxindex),) -> minindex == maxindex, indices)
            return indices
        end
    end

    return indices
end

function certifyindices(
    domain::SphericalTriangle{T},
    λs::AbstractVector{arb},
    vertex;
    show_trace = false,
) where {T}
    domain = SphericalTriangle(
        domain.angles[[mod1(i, 3) for i = vertex:vertex+2]],
        domain.parent,
    )

    # Define interval to work on and the parameters to use
    start = domain.parent(0)
    stop = -0.5 + sqrt(0.25 + getinterval(λs[end])[2] + 1)

    z = cos(theta_bound(domain, :upper))

    # Find all eigenvalues of the enclosing spherical cap sector in
    # the above interval
    eigenvalues = arb[]
    eigenvaluesflags = Int[]
    k = 1

    # Temporary for plotting
    νs = range(Float64(start), stop = Float64(stop), length = 2000)
    #λs = νs .* (νs .+ 1)
    res = []

    while true
        if T == fmpq
            μ = domain.parent(-k * inv(domain.angles[1]))
        else
            μ = -k * inv(domain.angles[1]) * domain.parent(π)
        end

        f = ν -> legendre_p(ν, μ, z)
        roots, flags = isolateroots(
            f,
            start,
            stop,
            evaltype = :taylor,
            maxevals = 5000,
            show_trace = false,
        )

        if !all(flags)
            @error "unable to isolate the roots of the legendre function"
        end

        push!(res, Float64.(f.(domain.parent.(νs))))

        if show_trace
            @show k
            @show length(roots)
        end

        if isempty(roots)
            break
        end

        append!(
            eigenvalues,
            [setinterval(l, u) * (setinterval(l, u) + 1) for (l, u) in roots],
        )
        append!(eigenvaluesflags, flags)
        k += 1
    end

    indices = Tuple{arb,Int,Int}[]
    for λ in λs
        minindex = length(filter(eig -> eig < λ, λs)) + 1

        # Find all eigenvalues of the enclosing domain that at least
        # as small as λ. Ensure that all of them have flag = 1
        idx = findall(eig -> !(eig > λ), eigenvalues)
        if !all(==(1), eigenvaluesflags[idx])
            @error "unable to upper bound the index for λ = $λ"
            maxindex = typemax(Int)
        else
            maxindex = length(idx)
        end
        push!(indices, (λ, minindex, maxindex))
    end

    return indices
    return [(λ, Int[length(filter(eig -> !(eig > λ), eigenvalues))]) for λ in λs]

    eigenvalues, λs, res
    # Compare with the given eigenvalues
end
