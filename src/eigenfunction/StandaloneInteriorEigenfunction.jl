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
    ks::UnitRange{Int};
    boundary = nothing,
    notransform::Bool = false,
) where {T <: Union{arb, arb_series}}
    if !notransform
        xy = coordinate_transformation(u, xy)
    end

    r, θ = polar_from_cartesian(xy)
    νs = u.stride*(div(ks.start, 2):div(ks.stop, 2))
    bessel_js = let sqrtλr = sqrt(λ)*r
        Dict(ν => bessel_j(u.parent(ν), sqrtλr) for ν in νs)
    end

    res = similar(ks, T)
    for i in eachindex(ks)
        k = ks[i]
        ν = u.stride*div(k, 2)

        if k == 1
            res[i] = bessel_js[ν]
        elseif iseven(k)
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
