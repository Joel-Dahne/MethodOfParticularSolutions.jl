function Base.show(io::IO, u::LinkedEigenfunction)
    println(io, "$(length(u.us)) linked $(eltype(u.us))")
    if !haskey(io, :compact) || !io[:compact]
        recur_io = IOContext(io, :compact => true)
        for v in u.us
            print(io, "- ")
            show(recur_io, v)
        end
        print(io, "number of set coefficients: $(length(coefficients(u)))")
    end
end

coefficients(u::LinkedEigenfunction) = coefficients(first(u.us))

function set_domain!(u::LinkedEigenfunction, domain::AbstractDomain)
    foreach(v -> set_domain!(v, domain), u.us)
    return u
end

function recompute!(u::LinkedEigenfunction)
    foreach(recompute!, u.us)
    return u
end

function set_eigenfunction!(
    u::LinkedEigenfunction{<:AbstractPlanarEigenfunction},
    coefficients::Vector,
)
    foreach(v -> set_eigenfunction!(v, coefficients), u.us)
    return u
end

function (u::LinkedEigenfunction)(
    xy::AbstractVector{T},
    λ::arb,
    ks::UnitRange{Int};
    boundary = nothing,
    notransform::Bool = false,
) where {T<:Union{arb,arb_series}}
    res = [zero(first(xy)) for _ in eachindex(ks)]

    for i in eachindex(u.us)
        if !(boundary ∈ u.excluded_boundaries[i])
            res .+=
                u.extra_coefficients[i] .*
                u.us[i](xy, λ, ks, boundary = boundary, notransform = notransform)
        end
    end

    return res
end
