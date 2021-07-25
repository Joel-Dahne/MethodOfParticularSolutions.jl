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

function set_eigenfunction!(u::LinkedEigenfunction, coefficients::Vector)
    foreach(v -> set_eigenfunction!(v, coefficients), u.us)
    return u
end

function (u::LinkedEigenfunction{S,T})(
    xy::AbstractVector,
    λ::Union{Real,arb},
    ks::UnitRange{Int};
    boundary = nothing,
    notransform::Bool = false,
) where {S,T}
    # It would be convenient to rely on the eigenfunctions in u.us to
    # convert the input to the appropriate type. However this is
    # problematic if all of them are excluded on the current boundary.
    # For performance reasons it's also better to do the conversion
    # only once.
    if S == arb
        if eltype(xy) == arb_series
            xy = convert(SVector{2,arb_series}, xy)
        else
            xy = convert(SVector{2,arb}, first(u.us).parent.(xy))
        end
        λ = first(u.us).parent(λ)
    else
        xy = convert(SVector{2,S}, xy)
        λ = convert(S, λ)
    end

    res = [zero(first(xy)) for _ in eachindex(ks)]

    for i in eachindex(u.us)
        if !(boundary ∈ u.excluded_boundaries[i])
            res += u.extra_coefficients[i] .* u.us[i](xy, λ, ks; boundary, notransform)
        end
    end

    return res
end
