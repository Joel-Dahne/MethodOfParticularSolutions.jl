StandaloneLightningEigenfunction(
    domain::AbstractPlanarDomain{S,T},
    i::Integer;
    l = 1,
    σ = 4,
    even::Bool = false,
    odd::Bool = false,
    reversed::Bool = false,
    outside = false,
) where {S,T} =
    StandaloneLightningEigenfunction{S,T}(domain, i; outside, l, σ, even, odd, reversed)

function StandaloneLightningEigenfunction{S,T}(
    domain::AbstractPlanarDomain,
    i::Integer;
    l = 1,
    σ = 4,
    even::Bool = false,
    odd::Bool = false,
    reversed::Bool = false,
    outside::Bool = false,
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

    return StandaloneLightningEigenfunction{S,T}(
        vertex(domain, i),
        orient,
        θ;
        l,
        σ,
        even,
        odd,
        reversed,
    )
end

function Base.show(io::IO, u::StandaloneLightningEigenfunction{T}) where {T}
    println(
        io,
        "StandaloneLightningEigenfunction{$T}" *
        ifelse(u.even, " - even", "") *
        ifelse(u.odd, " - odd", "") *
        ifelse(u.reversed, " - reversed", ""),
    )
    if !haskey(io, :compact) || !io[:compact]
        println(io, "vertex: $(u.vertex)")
        println(io, "orientation: $(u.orientation), θ: $(u.θ)")
        println(io, "l = $(u.l), σ = $(u.σ)")
        print(io, "number of set coefficients: $(length(u.coefficients))")
    end
end

"""
    coordinate_transformation(u::StandaloneLightningEigenfunction, xy::AbstractVector)

Takes a 2-element vector `xy` representing a point in the plane in
Cartesian coordinates and makes a (affine) change of coordinates so
that `u.vertex` is put at the origin and rotated `u.orientation + u.θ/2`
clockwise.
"""
function coordinate_transformation(
    u::StandaloneLightningEigenfunction{S,T},
    xy::AbstractVector,
) where {S,T}
    if T == fmpq
        s, c = convert.(S, sincospi(-u.orientation - u.θ // 2, u.parent))
    elseif T <: Rational
        s, c = sincospi(convert(T, -u.orientation - u.θ // 2))
    else
        s, c = sincos(-u.orientation - u.θ / 2)
    end

    M = SMatrix{2,2}(c, s, -s, c)
    res = M * (xy - u.vertex)

    if u.reversed
        return SVector(res[1], -res[2])
    else
        return res
    end
end

"""
    chargedistance(u::StandaloneLightningEigenfunction, i::Integer, n::Integer)

Compute the distance from the vertex for the `i`th charge point
assuming a total of `n` charge points are used.

We use a tapered distribution as proposed by Trefethen but we reverse
the order.

An old alternative version is `u.l*exp(-u.σ*(i - 1)/sqrt(u.parent(n)))`.
"""
chargedistance(u::StandaloneLightningEigenfunction{arb}, i::Integer, n::Integer) =
    u.l * exp(-u.σ * (sqrt(u.parent(n)) - sqrt(u.parent(n + 1 - i))))


chargedistance(u::StandaloneLightningEigenfunction{T}, i::Integer, n::Integer) where {T} =
    u.l * exp(-u.σ * (sqrt(convert(T, n)) - sqrt(convert(T, n + 1 - i))))

"""
    charge(
        u::StandaloneLightningEigenfunction,
        i::Integer,
        n::Integer,
        standard_coordinates = false,
)

Compute the location of the `i`th charge point assuming a total of `n`
points is used. The charge points are numbered from `1` to `n`.

They are returned in the coordinate system used by `u`, so they all
lie on the negative part of the x-axis. If `standard_coordinates` is
true then return them in the standard coordinate system instead.
"""
function charge(
    u::StandaloneLightningEigenfunction{S,T},
    i::Integer,
    n::Integer,
    standard_coordinates = false,
) where {S,T}
    d = chargedistance(u, i, n)
    if standard_coordinates
        if T == fmpq
            s, c = sincospi(-u.orientation - u.θ // 2, u.parent)
        elseif T <: Rational
            s, c = sincospi(convert(T, -u.orientation - u.θ // 2))
        else
            s, c = sincos(-u.orientation - u.θ / 2)
        end

        return u.vertex - d .* SVector{2,T}(c, s)
    else
        SVector{2}(-d, zero(d))
    end
end

function (u::StandaloneLightningEigenfunction{S,T})(
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

    res = similar(ks, eltype(xy))
    i = 1
    if u.even
        charge_indices = (div(ks.start - 1, 2)+1):(div(ks.stop - 1, 2)+1)
    elseif u.odd
        charge_indices = ks
    else
        charge_indices = (div(ks.start - 1, 3)+1):(div(ks.stop - 1, 3)+1)
    end
    for charge_index in charge_indices
        # Perform the change of coordinates corresponding to the charge
        # point
        xy_local = xy - charge(u, charge_index, charge_indices.stop)
        r, θ = polar_from_cartesian(xy_local)

        rsqrtλ = r * sqrt(λ)
        b = bessel_y(one(λ), rsqrtλ)
        s, c = sincos(θ)

        if u.even
            if charge_index == charge_indices.start == charge_indices.stop
                # Special case when first and last charge are the same
                if mod1(ks.start, 2) == 1
                    res[i] = bessel_y(zero(λ), rsqrtλ)
                    i += 1
                end
                if mod1(ks.stop, 2) == 2
                    res[i] = b * c
                    i += 1
                end
            elseif charge_index == charge_indices.start
                # Might not want all terms from the first charge
                if mod1(ks.start, 2) == 1
                    res[i] = bessel_y(zero(λ), rsqrtλ)
                    i += 1
                end
                res[i] = b * c
                i += 1
            elseif charge_index == charge_indices.stop
                # Might not want all terms from the last charge
                res[i] = bessel_y(zero(λ), rsqrtλ)
                i += 1
                if mod1(ks.stop, 2) == 2
                    res[i] = b * c
                    i += 1
                end
            else
                res[i] = bessel_y(zero(λ), rsqrtλ)
                res[i+1] = b * c
                i += 2
            end
        elseif u.odd
            # We only take one expansion from each charge in this
            # case. So it's a lot simpler.
            res[i] = b * s
            i += 1
        else
            if charge_index == charge_indices.start == charge_indices.stop
                # Special case when first and last charge are the same
                if mod1(ks.start, 3) == 1
                    res[i] = bessel_y(zero(λ), rsqrtλ)
                    i += 1
                end
                if mod1(ks.start, 3) <= 2 && mod1(ks.stop, 3) >= 2
                    res[i] = b * s
                    i += 1
                end
                if mod1(ks.stop, 3) == 3
                    res[i] = b * c
                    i += 1
                end
            elseif charge_index == charge_indices.start
                # Might not want all terms from the first charge
                if mod1(ks.start, 3) == 1
                    res[i] = bessel_y(zero(λ), rsqrtλ)
                    i += 1
                end
                if mod1(ks.start, 3) <= 2
                    res[i] = b * s
                    i += 1
                end
                res[i] = b * c
                i += 1
            elseif charge_index == charge_indices.stop
                # Might not want all terms from the last charge
                res[i] = bessel_y(zero(λ), rsqrtλ)
                i += 1
                if mod1(ks.stop, 3) >= 2
                    res[i] = b * s
                    i += 1
                end
                if mod1(ks.stop, 3) >= 3
                    res[i] = b * c
                    i += 1
                end
            else
                res[i] = bessel_y(zero(λ), rsqrtλ)
                res[i+1] = b * s
                res[i+2] = b * c
                i += 3
            end
        end
    end

    return res
end

# TODO: Figure out how to handle this
set_domain!(u::StandaloneLightningEigenfunction, ::AbstractDomain) = u

recompute!(u::StandaloneLightningEigenfunction) = u
