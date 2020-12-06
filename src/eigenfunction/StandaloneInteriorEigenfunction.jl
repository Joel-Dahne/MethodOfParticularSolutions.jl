function StandaloneInteriorEigenfunction(
    vertex::AbstractVector{arb},
    stride::Integer = 1,
    parent::ArbField = parent(first(vertex)),
)
    return StandaloneInteriorEigenfunction(
        vertex,
        stride,
        arb[],
        parent,
    )
end

StandaloneInteriorEigenfunction(domain::AbstractPlanarDomain, stride::Integer = 1) =
    StandaloneInteriorEigenfunction(center(domain), stride, arb[], domain.parent)

function Base.show(io::IO, u::StandaloneInteriorEigenfunction)
    println(io, "Standalone interior eigenfunction")
    if !haskey(io, :compact) || !io[:compact]
        println(io, "vertex: $(u.vertex)")
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
function coordinate_transformation(u::StandaloneInteriorEigenfunction, xy::AbstractVector)
    return xy .- u.vertex
end

function (u::StandaloneInteriorEigenfunction)(
    xy::AbstractVector{T},
    λ::arb,
    k::Integer;
    boundary = nothing,
    notransform::Bool = false,
) where {T <: Union{arb, arb_series}}
    if !notransform
        xy = coordinate_transformation(u, xy)
    end

    k = 1 + (k - 1)*u.stride

    r = sqrt(xy[1]^2 + xy[2]^2)
    θ = atan(xy[2], xy[1])
    ν = u.parent(div(k, 2))

    if k == 1
        return bessel_j(ν, sqrt(λ)*r)
    elseif k % 2 == 0
        return bessel_j(ν, sqrt(λ)*r)*sin(ν*θ)
    else
        return bessel_j(ν, sqrt(λ)*r)*cos(ν*θ)
    end
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

    νs = [div(1 + (k - 1)*u.stride, 2) for k in ks]
    bessel_js = let sqrtλr = sqrt(λ)*r
        Dict(ν => bessel_j(u.parent(ν), sqrtλr) for ν in unique(νs))
    end

    res = similar(ks, T)
    for i in eachindex(ks)
        k = 1 + (ks[i] - 1)*u.stride
        ν = νs[i]

        if k == 1
            res[i] = bessel_js[ν]
        elseif k % 2 == 0
            res[i] = bessel_js[ν]*sin(ν*θ)
        else
            res[i] = bessel_js[ν]*cos(ν*θ)
        end
    end

    return res
end

# TODO: Figure out how to handle this
set_domain!(u::StandaloneInteriorEigenfunction, ::AbstractDomain) = u

# TODO: Figure out how to handle this
recompute!(u::StandaloneInteriorEigenfunction) = u
