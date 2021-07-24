StandaloneVertexEigenfunction(
    domain::AbstractPlanarDomain{S,T},
    i::Integer;
    stride::Integer = 1,
    offset::Integer = 0,
    reversed::Bool = false,
    outside::Bool = false,
) where {S,T} =
    StandaloneVertexEigenfunction{S,T}(domain, i; stride, offset, reversed, outside)

function StandaloneVertexEigenfunction{S,T}(
    domain::AbstractPlanarDomain,
    i::Integer;
    stride::Integer = 1,
    offset::Integer = 0,
    reversed::Bool = false,
    outside = false,
) where {S,T}
    if T <: Union{Rational,fmpq} && !has_rational_angles(domain)
        throw(
            ErrorException(
                "can't construct an eigenfunction with rational angles from a domain with non-rational angles",
            ),
        )
    end

    if T <: Union{AbstractFloat,arb} && has_rational_angles(domain)
        θ = convert(T, angle(domain, i))
        orient = convert(T, orientation(domain, i))
    else
        # convert(Rational{U}, x::fmpq) only works when U is BigInt
        if T <: Rational && angle_raw(domain, i) isa fmpq
            θ = convert(T, convert(Rational{BigInt}, angle_raw(domain, i)))
            orient = convert(T, convert(Rational{BigInt}, orientation_raw(domain, i)))
        else
            θ = convert(T, angle_raw(domain, i))
            orient = convert(T, orientation_raw(domain, i))
        end
    end

    if outside
        orient += θ
        if has_rational_angles(domain)
            θ = 2 - θ
        else
            if S == arb
                θ = 2domain.parent(π) - θ
            else
                θ = 2convert(T, π) - θ
            end
        end
    end

    return StandaloneVertexEigenfunction{S,T}(
        vertex(domain, i),
        orient,
        θ;
        stride,
        offset,
        reversed,
        domain.parent,
    )
end

function Base.show(io::IO, u::StandaloneVertexEigenfunction)
    println(io, "Standalone vertex eigenfunction")
    if !haskey(io, :compact) || !io[:compact]
        println(io, "vertex: $(u.vertex)")
        println(io, "orientation: $(u.orientation)")
        println(io, "θ: $(u.θ)")
        print(io, "number of set coefficients: $(length(u.coefficients))")
    end
end

function set_eigenfunction!(u::StandaloneVertexEigenfunction, coefficients::Vector)
    resize!(u.coefficients, length(coefficients))
    copy!(u.coefficients, u.parent.(coefficients))
end

"""
    nu(u::StandaloneVertexEigenfunction, k::Integer)

Return `k*θ` as an `arb`, the parameter used for the Bessel function.
"""
nu(u::StandaloneVertexEigenfunction{T,fmpq}, k::Integer = 1) where {T} =
    u.parent(k * inv(u.θ))
nu(u::StandaloneVertexEigenfunction{T,arb}, k::Integer = 1) where {T} =
    u.parent(k * u.parent(π) / u.θ)

"""
    coordinate_transformation(u::StandaloneVertexEigenfunction, xy::AbstractVector)

Takes a 2-element vector `xy` representing a point in the plane in
Cartesian coordinates and makes a (affine) change of coordinates so
that `u.vertex` is put at the origin and rotated `u.orientation`
clockwise.
"""
function coordinate_transformation(
    u::StandaloneVertexEigenfunction{T,S},
    xy::AbstractVector,
) where {T,S}
    if S == fmpq
        if u.reversed
            s, c = sincospi(-u.orientation - u.θ, u.parent)
        else
            s, c = sincospi(-u.orientation, u.parent)
        end
    elseif S == arb
        s, c = sincos(-u.orientation)
    end
    M = SMatrix{2,2}(c, s, -s, c)
    res = M * (xy .- u.vertex)
    if u.reversed
        return SVector(res[1], -res[2])
    else
        return res
    end
end

function (u::StandaloneVertexEigenfunction)(
    xy::AbstractVector{T},
    λ::arb,
    ks::UnitRange{Int};
    boundary = nothing,
    notransform::Bool = false,
) where {T<:Union{arb,arb_series}}
    if !notransform
        xy = coordinate_transformation(u, xy)
    end

    r, θ = polar_from_cartesian(xy)

    # TODO: We have to choose a branch to work on. This does depend on
    # the domain but for now we only implement the one with θ on the
    # interval [0, 2π). Nemo doesn't implement mod2pi so we just do a
    # partial solution of adding 2π if it's below 0.
    if (T == arb && isnegative(θ)) || (T == arb_series && isnegative(θ[0]))
        if T == arb
            θ += 2parent(θ)(π)
        else
            θ += 2base_ring(parent(θ.poly))(π)
        end
    end

    rsqrtλ = r * sqrt(λ)
    res = similar(ks, T)
    for i in eachindex(ks)
        k = 1 + (ks[i] - 1) * u.stride + u.offset
        ν = nu(u, k)

        res[i] = bessel_j(ν, rsqrtλ) * sin(ν * θ)
    end

    return res
end

# TODO: Figure out how to handle this
set_domain!(u::StandaloneVertexEigenfunction, ::AbstractDomain) = u

# TODO: Figure out how to handle this
recompute!(u::StandaloneVertexEigenfunction) = u

# Have a function that computes u from the domain?
