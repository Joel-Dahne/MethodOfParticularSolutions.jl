function StandaloneInteriorEigenfunction(
    vertex::AbstractVector{arb},
    orientation::T = fmpq(0),
    parent::ArbField = parent(first(vertex));
    stride::Integer = 1,
    even = false,
) where {T <: Union{Rational,fmpq,arb}}
    if T == fmpq || T <: Rational
        orientation = fmpq(orientation)
    end
    return StandaloneInteriorEigenfunction(
        vertex,
        orientation,
        stride,
        even,
        arb[],
        parent,
    )
end

StandaloneInteriorEigenfunction(
    domain::AbstractPlanarDomain,
    orientation = fmpq(0);
    stride::Integer = 1,
    even = false,
) = StandaloneInteriorEigenfunction(center(domain), orientation, domain.parent; stride, even)

function Base.show(io::IO, u::StandaloneInteriorEigenfunction)
    println(
        io,
        "Standalone interior eigenfunction" * ifelse(u.even, " - even", ""),
    )
    if !haskey(io, :compact) || !io[:compact]
        println(io, "vertex: $(u.vertex)")
        println(io, "orientation: $(u.orientation)")
        print(io, "number of set coefficients: $(length(u.coefficients))")
    end
end

function set_eigenfunction!(u::StandaloneInteriorEigenfunction, coefficients::Vector)
    resize!(u.coefficients, length(coefficients))
    copy!(u.coefficients, u.parent.(coefficients))
end

"""
    coordinate_transformation(u::StandaloneInteriorEigenfunction, xy::AbstractVector)

Takes a 2-element vector `xy` representing a point in the plane in
Cartesian coordinates and makes a (affine) change of coordinates so
that `u.vertex` is put at the origin.
"""
function coordinate_transformation(
    u::StandaloneInteriorEigenfunction{T},
    xy::AbstractVector,
) where {T}
    iszero(u.orientation) && return xy .- u.vertex

    if T == fmpq
        s, c = sincospi(u.orientation, u.parent)
    elseif T == arb
        s, c = sincos(u.orientation)
    end
    M = SMatrix{2, 2}(c, s, -s, c)
    return M*(xy .- u.vertex)
end

function (u::StandaloneInteriorEigenfunction)(
    xy::AbstractVector{T},
    λ::arb,
    ks::UnitRange{Int};
    boundary = nothing,
    notransform::Bool = false,
) where {T <: Union{arb, arb_series}}
    if !notransform
        xy = coordinate_transformation(u, xy)
    end

    r, θ = polar_from_cartesian(xy)
    if u.even
        νs = u.stride*(ks.start - 1:ks.stop - 1)
    else
        νs = u.stride*(div(ks.start, 2):div(ks.stop, 2))
    end
    bessel_js = let sqrtλr = sqrt(λ)*r
        Dict(ν => bessel_j(u.parent(ν), sqrtλr) for ν in νs)
    end

    res = similar(ks, T)
    for i in eachindex(ks)
        k = ks[i]
        if u.even
            ν = u.stride*(k - 1)

            res[i] = bessel_js[ν]*cos(ν*θ)
        else
            ν = u.stride*div(k, 2)

            if k == 1
                res[i] = bessel_js[ν]
            elseif iseven(k)
                res[i] = bessel_js[ν]*sin(ν*θ)
            else
                res[i] = bessel_js[ν]*cos(ν*θ)
            end
        end
    end

    return res
end

# TODO: Figure out how to handle this
set_domain!(u::StandaloneInteriorEigenfunction, ::AbstractDomain) = u

# TODO: Figure out how to handle this
recompute!(u::StandaloneInteriorEigenfunction) = u
