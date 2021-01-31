export PayneEigenfunction

struct PayneEigenfunction <: AbstractPlanarEigenfunction
    u::CombinedEigenfunction
    weights::Vector{Int}

    """
        PayneEigenfunction(u::CombinedEigenfunction; weights = [1, 1, 1, 4])

    This is a type mainly used for dispatch. We do some specific
    things for the Payne polygon and this is captured by this type. It
    does some simple checks to make sure that the given eigenfunction
    has the assumed properties.

    `weights` determines the distribution of collocation points
    between the different expansions. `weights = [1, 1, 1, 1]` would
    mean that the collocation points are spread out between the
    expansions only depending on how many free coefficients each one
    has. `weights = [1, 1, 1, 4]` means that the fourth expansion gets
    four times as many collocation points per free coefficient.
    """
    function PayneEigenfunction(
        u::CombinedEigenfunction;
        weights = [1, 1, 1, 4],
    )
        @assert length(u.us) == 4
        @assert length(weights) == 4

        @assert u.us[1] isa LinkedEigenfunction

        @assert u.us[2] isa LinkedEigenfunction{<:StandaloneLightningEigenfunction}
        @assert u.us[2].us[1].even
        @assert !u.us[2].us[1].reversed

        @assert u.us[3] isa LinkedEigenfunction{<:StandaloneLightningEigenfunction}
        @assert !u.us[3].us[1].even
        @assert !u.us[3].us[7].even
        @assert !u.us[3].us[1].reversed
        @assert u.us[3].us[7].reversed

        @assert u.us[4] isa StandaloneInteriorEigenfunction{fmpq}
        @assert u.us[4].stride == 6
        @assert u.us[4].even

        return new(u, weights)
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

            pts = map(d -> vertex(domain, 1) + d.*w, distances)
            return pts, fill(1, length(pts))
        elseif i == 2
            w = vertex(domain, 9) - vertex(domain, 8)
            @assert !(δs[1] > LinearAlgebra.norm(w))
            w = normalize(w)

            distances = vcat(
                [[j*δ/M for j in 1:M] for (δ, M) in zip(δs, per_charge)]...
            )
            distances = sort(distances, by = Float64, rev = true)

            pts = map(d -> vertex(domain, 8) + d.*w, distances)
            return pts, fill(8, length(pts))
        elseif i == 3
            # Distribute the points between the two sides
            left_side = cld.(per_charge, 2)
            right_side = div.(per_charge, 2)

            w₁ = vertex(domain, 8) - vertex(domain, 9)
            w₂ = vertex(domain, 7) - vertex(domain, 9)
            @assert !(δs[1] > LinearAlgebra.norm(w₁))
            @assert !(δs[1] > LinearAlgebra.norm(w₂))
            w₁ = normalize(w₁)
            w₂ = normalize(w₂)

            distances₁ = vcat(
                [[j*δ/M for j in 1:M] for (δ, M) in zip(δs, left_side)]...
            )
            distances₁ = sort(distances₁, by = Float64, rev = true)
            distances₂ = vcat(
                [[j*δ/M for j in 1:M] for (δ, M) in zip(δs, right_side)]...
            )
            distances₂ = sort(distances₂, by = Float64, rev = true)

            pts₁ = map(d -> vertex(domain, 9) + d.*w₁, distances₁)
            pts₂ = map(d -> vertex(domain, 9) + d.*w₂, distances₂)
            return vcat(pts₁, pts₂), vcat(fill(8, length(pts₁)), fill(9, length(pts₂)))
        end
    elseif i == 4
        return boundary_points(domain, u, N, n)
    else
        throw(ArgumentError("invalid value i = $i"))
    end
end

function boundary_points(
    domain::AbstractDomain,
    u::PayneEigenfunction,
    N::Integer,
    n::Integer;
    distribution = :mixed,
)
    if distribution == :chebyshev
        # Distribute points normally
        return boundary_points(domain, u.u, N, n; distribution)
    elseif distribution == :mixed
        # Distribute the n collocations points between the different
        # expansions
        ns = let
            # Compute weights for the different expansions
            ws = u.u.orders .* u.weights

            # Distribute collocation points according to the weights We do
            # it "backwards" to prioritize the interior expansion
            fullcycles, R = divrem(n, sum(ws))
            A = reverse([0; cumsum(reverse(ws))])[2:end]
            remaining = [max(min(ws[i], R - A[i]), 0) for i in eachindex(ws)]
            fullcycles*ws + remaining
        end

        return boundary_points(domain, u, N, tuple(ns...); distribution)
    else
        throw(ArgumentError("invalid distribution $distribution"))
    end
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
    res = vcat(res...)
    @assert sum(length, r[1] for r in res) == sum(ns)
    return vcat(getindex.(res, 1)...), vcat(getindex.(res, 2)...)
end
