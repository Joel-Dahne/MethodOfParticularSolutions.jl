StandaloneInteriorEigenfunction(
    domain::AbstractPlanarDomain{S,T},
    orientation = T == arb ? domain.parent(0) : zero(T);
    stride::Integer = 1,
    offset::Integer = 0,
    even = false,
    odd = false,
) where {S,T} = StandaloneInteriorEigenfunction(
    center(domain),
    T == arb ? domain.parent(orientation) : convert(T, orientation);
    stride,
    offset,
    even,
    odd,
    domain.parent,
)

StandaloneInteriorEigenfunction{S,T}(
    domain::AbstractPlanarDomain,
    orientation = zero(T);
    stride::Integer = 1,
    offset::Integer = 0,
    even = false,
    odd = false,
) where {S,T} = StandaloneInteriorEigenfunction{S,T}(
    center(domain),
    convert(T, orientation);
    stride,
    offset,
    even,
    odd,
    domain.parent,
)

function Base.show(io::IO, u::StandaloneInteriorEigenfunction)
    println(
        io,
        "Standalone interior eigenfunction" *
        ifelse(u.even, " - even", "") *
        ifelse(u.odd, " - odd", ""),
    )
    if !haskey(io, :compact) || !io[:compact]
        println(io, "vertex: $(u.vertex)")
        println(io, "orientation: $(u.orientation)")
        print(io, "number of set coefficients: $(length(u.coefficients))")
    end
end

"""
    coordinate_transformation(u::StandaloneInteriorEigenfunction, xy::AbstractVector)

Takes a 2-element vector `xy` representing a point in the plane in
Cartesian coordinates and makes a (affine) change of coordinates so
that `u.vertex` is put at the origin and rotated `u.orientation`
clockwise.
"""
function coordinate_transformation(
    u::StandaloneInteriorEigenfunction{S,T},
    xy::AbstractVector,
) where {S,T}
    iszero(u.orientation) && return xy - u.vertex

    if T == fmpq
        s, c = convert.(S, sincospi(-u.orientation, u.parent))
    elseif T <: Rational
        s, c = sincospi(convert(S, -u.orientation))
    else
        s, c = sincos(-u.orientation)
    end

    M = SMatrix{2,2}(c, s, -s, c)
    return M * (xy - u.vertex)
end

function (u::StandaloneInteriorEigenfunction{S,T})(
    xy::AbstractVector,
    λ::Union{Real,arb},
    ks::UnitRange{Int};
    boundary = nothing,
    notransform::Bool = false,
) where {S,T}
    # Convert input to type S. One exception is of xy is of type arb_series
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

    if !notransform
        xy = coordinate_transformation(u, xy)
    end

    r, θ = polar_from_cartesian(xy)

    if u.even
        νs = u.stride * (ks.start-1:ks.stop-1) .+ u.offset
    elseif u.odd
        νs = u.stride * (ks.start-1:ks.stop-1) .+ 1 .+ u.offset
    else
        νs = u.stride * (div(ks.start, 2):div(ks.stop, 2)) .+ u.offset
    end
    bessel_js = let sqrtλr = sqrt(λ) * r
        if S == arb
            Dict(ν => bessel_j(u.parent(ν), sqrtλr) for ν in νs)
        else
            Dict(ν => bessel_j(convert(S, ν), sqrtλr) for ν in νs)
        end
    end

    res = similar(ks, eltype(xy))
    for i in eachindex(ks)
        k = ks[i]
        if u.even
            ν = u.stride * (k - 1) + u.offset

            res[i] = bessel_js[ν] * cos(ν * θ)
        elseif u.odd
            ν = u.stride * (k - 1) + 1 + u.offset

            res[i] = bessel_js[ν] * sin(ν * θ)
        else
            ν = u.stride * div(k, 2) + u.offset

            if isodd(k)
                res[i] = bessel_js[ν] * cos(ν * θ)
            else
                res[i] = bessel_js[ν] * sin(ν * θ)
            end
        end
    end

    return res
end
