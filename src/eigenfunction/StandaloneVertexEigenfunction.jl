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
    nu(u::StandaloneVertexEigenfunction{S,T}, k::Integer)

Return `k / (angle(u, k) / π)` converted to type `S`. This is the
parameter used for the Bessel function.
"""
nu(
    u::StandaloneVertexEigenfunction{S,<:Rational},
    k::Integer = 1,
) where {S<:AbstractFloat} = convert(S, k / u.θ)
# Cant convert fmpq to S in general
nu(u::StandaloneVertexEigenfunction{S,fmpq}, k::Integer = 1) where {S<:AbstractFloat} =
    convert(S, convert(Rational{BigInt}, k * inv(u.θ)))
nu(u::StandaloneVertexEigenfunction{S,T}, k::Integer = 1) where {S<:AbstractFloat,T} =
    convert(S, k * (π / u.θ))
# Use parent for converting to arb
nu(u::StandaloneVertexEigenfunction{arb,<:Union{Rational,fmpq}}, k::Integer = 1) =
    u.parent(k * inv(u.θ))
nu(u::StandaloneVertexEigenfunction{arb,T}, k::Integer = 1) where {T} =
    u.parent(k * u.parent(π) / u.θ)

"""
    coordinate_transformation(u::StandaloneVertexEigenfunction, xy::AbstractVector)

Takes a 2-element vector `xy` representing a point in the plane in
Cartesian coordinates and makes a (affine) change of coordinates so
that `u.vertex` is put at the origin and rotated `u.orientation`
clockwise.

If `u.reversed = true` then it is instead rotated `u.orientation +
u.θ` and also flips the sign of the `y`-axis.
"""
function coordinate_transformation(
    u::StandaloneVertexEigenfunction{S,T},
    xy::AbstractVector,
) where {S,T}
    angle = u.reversed ? u.orientation + u.θ : u.orientation

    if T == fmpq
        s, c = sincospi(-angle, u.parent)
    elseif T <: Rational
        s, c = sincospi(convert(S, -angle))
    else
        s, c = sincos(-angle)
    end

    M = SMatrix{2,2}(c, s, -s, c)
    res = M * (xy - u.vertex)

    if u.reversed
        return SVector(res[1], -res[2])
    else
        return res
    end
end

function (u::StandaloneVertexEigenfunction{S,T})(
    xy::AbstractVector,
    λ::Union{Real,arb},
    ks::UnitRange{Int};
    boundary = nothing,
    notransform::Bool = false,
) where {S,T}
    # Promote to common type. Neither arb nor arb_series supports
    # promote so these we handle separately.
    if S == arb
        if eltype(xy) == arb_series
            U = arb_series
            xy = convert(SVector{2,arb_series}, xy)
        else
            U = arb
            xy = convert(SVector{2,arb}, u.parent.(xy))
        end
        λ = u.parent(λ)
    else
        U = promote_type(S, eltype(xy), typeof(λ))
        xy = convert(SVector{2,U}, xy)
        # TODO: Here we will eventually have to handle ArbSeries
        # differently
        λ = convert(U, λ)
    end

    if !notransform
        xy = coordinate_transformation(u, xy)
    end

    r, θ = polar_from_cartesian(xy)

    # TODO: We have to choose a branch to work on in θ. This does
    # depend on the domain but for now we only implement the one with
    # θ on the interval [0, 2π). Nemo doesn't implement mod2pi so we
    # just do a partial solution of adding 2π if it's below 0.
    if (U == arb_series && θ[0] < 0) || (U != arb_series && θ < 0)
        if U == arb_series
            θ += 2base_ring(parent(θ.poly))(π)
        elseif U == arb
            θ += 2parent(θ)(π)
        else
            θ += 2convert(U, π)
        end
    end

    rsqrtλ = r * sqrt(λ)
    res = similar(ks, U)
    for i in eachindex(ks)
        k = 1 + (ks[i] - 1) * u.stride + u.offset
        ν = nu(u, k)

        res[i] = bessel_j(ν, rsqrtλ) * sin(ν * θ)
    end

    return res
end

# TODO: Figure out how to handle this
set_domain!(u::StandaloneVertexEigenfunction, ::AbstractDomain) = u

recompute!(u::StandaloneVertexEigenfunction) = u
