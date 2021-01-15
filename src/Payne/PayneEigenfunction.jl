export PayneEigenfunction

struct PayneEigenfunction <: AbstractPlanarEigenfunction
    u::CombinedEigenfunction

    """
        PayneEigenfunction(u::CombinedEigenfunction)

    This is a type mainly used for dispatch. We do some specific
    things for the Payne polygon and this is captured by this type. It
    does some simple checks to make sure that the given eigenfunction
    has the assumed properties.
    """
    function PayneEigenfunction(u::CombinedEigenfunction)
        @assert u.us[1].us[1].even
        @assert !u.us[2].us[1].even && u.us[2].us[7].reversed
        @assert !u.us[3].us[1].even && u.us[3].us[7].reversed
        @assert ((u.us[4] isa LinkedEigenfunction) && u.us[4].us[1].even) ||
            ((u.us[4] isa StandaloneInteriorEigenfunction) && u.us[4].even)
        return new(u)
    end

end

coefficients(u::PayneEigenfunction) = coefficients(u.u)

function set_domain!(u::PayneEigenfunction, domain::AbstractDomain)
    set_domain!(u.u, domain)
    return u
end

function recompute!(u::PayneEigenfunction)
    recompute!(u.u)
    return u
end

function set_eigenfunction!(
    u::PayneEigenfunction,
    coefficients::Vector,
)
    set_eigenfunction!(u.u, coefficients)
    return u
end

function (u::PayneEigenfunction)(
    xy::AbstractVector{T},
    λ::arb,
    ks::UnitRange{Int};
    boundary = nothing,
    notransform::Bool = false,
) where {T<:Union{arb,arb_series}}
    return u.u(xy, λ, ks; boundary, notransform)
end

"""
    payne_boundary(domain, u::PayneEigenfunction, N, i, n)

Return `n` collocation points according to the `i`th expansion.
"""
function payne_boundary(domain, u::PayneEigenfunction, N, i, n)
    iszero(n) && return SVector{2,arb}[], Int[]
    u = u.u
    if i == 1 || i == 2 || i == 3
        # Compute number of charges used
        # FIXME: This is not always correct
        num_charges = cld(u.orders[i]*N, ifelse(i == 1, 2, 3)*sum(u.orders))
        δs = [chargedistance(u.us[i].us[1], j, n) for j in 1:num_charges]

        # Distribute the points over the charges
        per_charge = fill(div(n, num_charges), num_charges) +
            [ifelse(n%num_charges >= j, 1, 0) for j in 1:num_charges]

        if i == 1
            w = vertex(domain, 2) - vertex(domain, 1)
            @assert !(δs[1] > LinearAlgebra.norm(w))
            w = normalize(w)

            distances = vcat(
                [[j*δ/M for j in 1:M] for (δ, M) in zip(δs, per_charge)]...
            )
            distances = sort(distances, by = Float64, rev = true)
            #@show Float64.(distances)
            pts = map(d -> vertex(domain, 1) + d.*w, distances)
            return pts, fill(1, length(pts))
        elseif i == 2
            # Distribute the points between the two sides
            left_side = cld.(per_charge, 2)
            right_side = div.(per_charge, 2)
            #@show left_side right_side

            w₁ = vertex(domain, 7) - vertex(domain, 8)
            w₂ = vertex(domain, 9) - vertex(domain, 8)
            @assert !(δs[1] > LinearAlgebra.norm(w₁))
            @assert !(δs[1] > LinearAlgebra.norm(w₂))
            w₁ = normalize(w₁)
            w₂ = normalize(w₂)

            distances₁ = vcat(
                [[j*δ/M for j in 1:M] for (δ, M) in zip(δs, left_side)]...
            )
            distances₁ = sort(distances₁, by = Float64, rev = true)
            #@show Float64.(distances₁)
            pts₁ = map(d -> vertex(domain, 8) + d.*w₁, distances₁)
            distances₂ = vcat(
                [[j*δ/M for j in 1:M] for (δ, M) in zip(δs, right_side)]...
            )
            distances₂ = sort(distances₂, by = Float64, rev = true)
            #@show Float64.(distances₂)
            pts₂ = map(d -> vertex(domain, 8) + d.*w₂, distances₂)
            return vcat(pts₁, pts₂), vcat(fill(7, length(pts₁)), fill(8, length(pts₂)))
        elseif i == 3
            # Distribute the points between the two sides
            left_side = cld.(per_charge, 2)
            right_side = div.(per_charge, 2)
            #@show left_side right_side

            w₁ = vertex(domain, 10) - vertex(domain, 7)
            w₂ = vertex(domain, 8) - vertex(domain, 7)
            @assert !(δs[1] > LinearAlgebra.norm(w₁))
            #@assert !(δs[1] > LinearAlgebra.norm(w₂))
            w₁ = normalize(w₁)
            w₂ = normalize(w₂)

            distances₁ = vcat(
                [[j*δ/M for j in 1:M] for (δ, M) in zip(δs, left_side)]...
            )
            distances₁ = sort(distances₁, by = Float64, rev = true)
            #@show Float64.(distances₁)
            pts₁ = map(d -> vertex(domain, 7) + d.*w₁, distances₁)
            distances₂ = vcat(
                [[j*δ/M for j in 1:M] for (δ, M) in zip(δs, right_side)]...
            )
            distances₂ = sort(distances₂, by = Float64, rev = true)
            #@show Float64.(distances₂)
            pts₂ = map(d -> vertex(domain, 7) + d.*w₂, distances₂)
            return vcat(pts₁, pts₂), vcat(fill(10, length(pts₁)), fill(7, length(pts₂)))
        end
    elseif i == 4
        # Distribute the points over the boundaries
        per_boundary = fill(div(n, 4), 4) + [ifelse(n%4 >= j, 1, 0) for j in 1:4]
        res = [
            boundary_points(domain, 7, per_boundary[1]);
            [
                getindex.(boundary_points(domain, k, 2per_boundary[j + 1]), Ref(1:per_boundary[j + 1]))
                for (j, k) in enumerate([1, 8, 10])
            ]
        ]
        return vcat(getindex.(res, 1)...), vcat(getindex.(res, 2)...)
    else
        throw(ArgumentError("invalid value i = $i"))
    end
end

function boundary_points(
    domain::AbstractDomain,
    u::PayneEigenfunction,
    N::Integer,
    n::Integer;
    distribution = :chebyshev,
)
    # Distribute the n collocations points between the different
    # expansions
    ns = let u = u.u
        # Compute weights for the different expansions, this is just
        # the distribution of the free coefficients and giving some
        # extra weight to the interior expansion
        ws = u.orders .* [1, 1, 1, 4]

        # Distribute collocation points according to the weights We do
        # it "backwards" to prioritize the interior expansion
        A = [0; cumsum(reverse(ws))]
        fullcycles, R = divrem(n, A[end])
        remaining = reverse([max(min(ws[i], R - A[i]), 0) for i in reverse(eachindex(ws))])
        [fullcycles*ws[i] + remaining[i] for i in eachindex(ws)]
    end

    return boundary_points(domain, u, N, tuple(ns...); distribution)

    return let u = u.u
        m = 4
        per_side = fill(div(n, m), m) + [ifelse(n%m >= j, 1, 0) for j in 1:m]

        res = [
            boundary_points(domain, 7, per_side[1]; distribution);
            [
                getindex.(boundary_points(domain, i, 2M; distribution), Ref(1:M))
                for (i, M) in zip([1, 8, 10], per_side[2:end])
            ]
        ]

        vcat(getindex.(res, 1)...), vcat(getindex.(res, 2)...)
    end

    return boundary_points(domain, u.u, n; distribution)


end

"""
    boundary_points(domain, u::PaynePolygon, N::Integer, ns)

`n₁, n₂, n₃` should be the number of charges for these expansions.
`n₄` should be the number of free coefficients for the interior
expansion.

Take `ns[i]` collocation points according to the `i`th expansion.
"""
function boundary_points(
    domain::AbstractDomain,
    u::PayneEigenfunction,
    N::Integer,
    ns::NTuple{4,<:Integer};
    distribution = :chebyshev,
)
    res = [payne_boundary(domain, u, N, i, ns[i]) for i in 1:4]
    @show length(res[1][1]) length(res[2][1]) length(res[3][1]) length(res[4][1])
    res = vcat(res...)
    @assert sum(length, r[1] for r in res) == sum(ns)
    return vcat(getindex.(res, 1)...), vcat(getindex.(res, 2)...)
end
