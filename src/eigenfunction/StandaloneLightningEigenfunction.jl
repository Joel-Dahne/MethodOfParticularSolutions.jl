function StandaloneLightningEigenfunction(
    domain::Triangle{T},
    i::Integer;
    outside = false,
    l::arb = domain.parent(1),
    σ::arb = domain.parent(4),
    even::Bool = false,
    reversed::Bool = false,
) where {T <: Union{arb,fmpq}}
    θ = domain.angles[i]

    if i == 1
        orientation = zero(θ)
    elseif i == 2
        orientation = ifelse(T == arb, domain.parent(π), 1) - θ
    elseif i == 3
        orientation = 2ifelse(T == arb, domain.parent(π), 1) -
            domain.angles[2] - θ
    end

    if outside
        orientation += θ
        θ = ifelse(T == arb, 2domain.parent(π), 2) - θ
    end

    return StandaloneLightningEigenfunction(
        vertex(domain, i),
        orientation,
        θ,
        domain.parent;
        l,
        σ,
        even,
        reversed,
    )
end

# TODO: Currently the orientation is computed from the placement of
# the vertices and will therefore always be an arb and never and fmpq.
# The value for θ is also converted to an arb for that reason.
function StandaloneLightningEigenfunction(
    domain::Polygon{T},
    i::Integer;
    outside = false,
    l::arb = domain.parent(1),
    σ::arb = domain.parent(4),
    even::Bool = false,
    reversed::Bool = false,
) where {T <: Union{arb,fmpq}}
    θ = angle(domain, i)

    v = vertex(domain, mod1(i + 1, length(vertices(domain)))) - vertex(domain, i)
    orientation = atan(v[2], v[1])

    if contains_zero(v[2])
        @warn "orientation could not be computed accurately, orientation = $orientation"
    end

    if outside
        orientation += θ
        θ = 2domain.parent(π) - θ
    end

    return StandaloneLightningEigenfunction(
        vertex(domain, i),
        orientation,
        θ,
        domain.parent;
        l,
        σ,
        even,
        reversed,
    )
end

function StandaloneLightningEigenfunction(
    domain::TransformedDomain,
    i::Integer;
    outside = false,
    l::arb = domain.parent(1),
    σ::arb = domain.parent(4),
    even::Bool = false,
    reversed::Bool = false,
)
    u = StandaloneLightningEigenfunction(domain.original, i; outside, l, σ, even, reversed)
    if typeof(u.orientation) == typeof(domain.rotation)
        orientation = u.orientation + domain.rotation
    else
        if u.orientation isa fmpq
            orientation = domain.parent(π)*u.orientation + domain.rotation
        else
            orientation = u.orientation + domain.parent(π)*domain.rotation
        end
    end

    u = StandaloneLightningEigenfunction(
        domain.map(u.vertex),
        orientation,
        u.θ,
        domain.parent,
        l = u.l,
        σ = u.σ,
        even = u.even,
        reversed = u.reversed,
    )

    return u
end

function Base.show(io::IO, u::StandaloneLightningEigenfunction)
    println(
        io,
        "Standalone lightning eigenfunction" *
        ifelse(u.even, " - even", "") *
        ifelse(u.reversed, " - reversed", "")
        ,
    )
    if !haskey(io, :compact) || !io[:compact]
        println(io, "vertex: $(u.vertex)")
        println(io, "orientation: $(u.orientation)")
        println(io, "θ: $(u.θ)")
        println(io, "l = $(u.l), σ = $(u.σ)")
        print(io, "number of set coefficients: $(length(u.coefficients))")
    end
end

function set_eigenfunction!(u::StandaloneLightningEigenfunction, coefficients::Vector)
    resize!(u.coefficients, length(coefficients))
    copy!(u.coefficients, u.parent.(coefficients))
end

"""
    coordinate_transformation(u::StandaloneLightningEigenfunction, xy::AbstractVector)

Takes a 2-element vector `xy` representing a point in the plane in
Cartesian coordinates and makes a (affine) change of coordinates so
that `u.vertex` is put at the origin and rotated `u.orientation + u.θ/2`
clockwise.
"""
function coordinate_transformation(
    u::StandaloneLightningEigenfunction{T,S},
    xy::AbstractVector,
) where {T,S <: Union{fmpq,arb}}
    if S == fmpq
        s, c = sincospi(-u.orientation - u.θ//2, u.parent)
    elseif S == arb
        s, c = sincos(-u.orientation - u.θ/2)
    end
    M = SMatrix{2, 2}(c, s, -s, c)
    res = M*(xy .- u.vertex)
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
"""
chargedistance(u::StandaloneLightningEigenfunction, i::Integer, n::Integer)=
    u.l*exp(-u.σ*(sqrt(u.parent(n)) - sqrt(u.parent(n + 1 - i))))
# An old alternative version is
#u.l*exp(-u.σ*(i - 1)/sqrt(u.parent(n)))

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
    u::StandaloneLightningEigenfunction{T,S},
    i::Integer,
    n::Integer,
    standard_coordinates = false,
) where {T,S<:Union{arb,fmpq}}
    d = chargedistance(u, i, n)
    if standard_coordinates
        if T == fmpq
            s, c = sincospi(u.orientation + u.θ//2, u.parent)
        else
            s, c = sincos(u.orientation + u.θ/2)
        end
        return u.vertex - d.*SVector{2,arb}(c, s)
    else
        SVector{2,arb}(-d, zero(d))
    end
end

function (u::StandaloneLightningEigenfunction)(
    xy::AbstractVector{T},
    λ::arb,
    ks::UnitRange{Int};
    boundary = nothing,
    notransform::Bool = false,
    n = div(length(ks) - 1, 3) + 1, # Total number of charge points
) where {T <: Union{arb, arb_series}}
    if !notransform
        xy = coordinate_transformation(u, xy)
    end

    res = similar(ks, T)
    i = 1
    if u.even
        charge_indices = (div(ks.start-1, 2)+1):(div(ks.stop-1, 2)+1)
    else
        charge_indices = (div(ks.start-1, 3)+1):(div(ks.stop-1, 3)+1)
    end
    for charge_index in charge_indices
        # Perform the change of coordinates corresponding to the charge
        # point
        xy_local = xy - charge(u, charge_index, n)
        r, θ = polar_from_cartesian(xy_local)

        rsqrtλ = r*sqrt(λ)
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
                    res[i] = b*c
                    i += 1
                end
            elseif charge_index == charge_indices.start
                # Might not want all terms from the first charge
                if mod1(ks.start, 2) == 1
                    res[i] = bessel_y(zero(λ), rsqrtλ)
                    i += 1
                end
                res[i] = b*c
                i += 1
            elseif charge_index == charge_indices.stop
                # Might not want all terms from the last charge
                res[i] = bessel_y(zero(λ), rsqrtλ)
                i += 1
                if mod1(ks.stop, 2) == 2
                    res[i] = b*c
                    i += 1
                end
            else
                res[i] = bessel_y(zero(λ), rsqrtλ)
                res[i + 1] = b*c
                i += 2
            end
        else
            if charge_index == charge_indices.start == charge_indices.stop
                # Special case when first and last charge are the same
                if mod1(ks.start, 3) == 1
                    res[i] = bessel_y(zero(λ), rsqrtλ)
                    i += 1
                end
                if mod1(ks.start, 3) <= 2 && mod1(ks.stop, 3) >= 2
                    res[i] = b*s
                    i += 1
                end
                if mod1(ks.stop, 3) == 3
                    res[i] = b*c
                    i += 1
                end
            elseif charge_index == charge_indices.start
                # Might not want all terms from the first charge
                if mod1(ks.start, 3) == 1
                    res[i] = bessel_y(zero(λ), rsqrtλ)
                    i += 1
                end
                if mod1(ks.start, 3) <= 2
                    res[i] = b*s
                    i += 1
                end
                res[i] = b*c
                i += 1
            elseif charge_index == charge_indices.stop
                # Might not want all terms from the last charge
                res[i] = bessel_y(zero(λ), rsqrtλ)
                i += 1
                if mod1(ks.stop, 3) >= 2
                    res[i] = b*s
                    i += 1
                end
                if mod1(ks.stop, 3) >= 3
                    res[i] = b*c
                    i += 1
                end
            else
                res[i] = bessel_y(zero(λ), rsqrtλ)
                res[i + 1] = b*s
                res[i + 2] = b*c
                i += 3
            end
        end
    end

    return res
end

# TODO: Figure out how to handle this
set_domain!(u::StandaloneLightningEigenfunction, ::AbstractDomain) = u

# TODO: Figure out how to handle this
recompute!(u::StandaloneLightningEigenfunction) = u
