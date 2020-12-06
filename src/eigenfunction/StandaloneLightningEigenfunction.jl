function StandaloneLightningEigenfunction(
    vertex::SVector{2,arb},
    orientation::T,
    θ::T,
    parent::ArbField = parent(vertex[1]);
    l::arb = parent(1),
    σ::arb = parent(2.5),
) where {T <: Union{arb,fmpq}}
    return StandaloneLightningEigenfunction(
        vertex,
        orientation,
        θ,
        l,
        σ,
        arb[],
        parent,
    )
end

function StandaloneLightningEigenfunction(
    domain::Triangle{T},
    i::Integer;
    outside = false,
    l::arb = domain.parent(1),
    σ::arb = domain.parent(2.5),
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
        σ
    )
end

function StandaloneLightningEigenfunction(
    domain::Polygon{T},
    i::Integer;
    outside = false,
    l::arb = domain.parent(1),
    σ::arb = domain.parent(2.5),
) where {T <: Union{arb,fmpq}}
    θ = domain.angles[i]

    orientation = (i - 1)*ifelse(T == arb, domain.parent(π), 1) - θ

    if outside
        throw(ErrorException("Not implemented"))
    end

    return StandaloneLightningEigenfunction(
        vertex(domain, i),
        orientation,
        θ,
        domain.parent;
        l,
        σ
    )
end

function StandaloneLightningEigenfunction(
    domain::TransformedDomain,
    i::Integer;
    outside = false,
    l::arb = domain.parent(1),
    σ::arb = domain.parent(2.5),
)
    u = StandaloneLightningEigenfunction(domain.original, i; outside, l, σ)
    # FIXME: This only works if u.θ::arb or if u.θ and
    # domain.orientation both are fmpq.
    u = StandaloneLightningEigenfunction(
        domain.map(u.vertex),
        u.orientation + domain.rotation,
        u.θ,
        u.l,
        u.σ,
        u.coefficients,
        domain.parent
    )

    return u
end

function Base.show(io::IO, u::StandaloneLightningEigenfunction)
    println(io, "Standalone lightning eigenfunction")
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
function coordinate_transformation(u::StandaloneLightningEigenfunction{fmpq}, xy::AbstractVector)
    s, c = sincospi(-u.orientation - u.θ//2, u.parent)
    M = SMatrix{2, 2}(c, s, -s, c)
    return M*(xy .- u.vertex)
end

function coordinate_transformation(u::StandaloneLightningEigenfunction{arb}, xy::AbstractVector)
    s, c = sincos(-u.orientation - u.θ/2)
    M = SMatrix{2, 2}(c, s, -s, c)
    return M*(xy .- u.vertex)
end

"""
    chargedistance(u::StandaloneLightningEigenfunction, i::Integer, n::Integer)

Compute the distance from the vertex for the `i`th charge point
assuming a total of `n` charge points are used.
"""
chargedistance(u::StandaloneLightningEigenfunction, i::Integer, n::Integer)=
    u.l*exp(-u.σ*(i - 1)/sqrt(u.parent(n)))

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
lie on the negative part of the x-axis. If `defualt_coordinates` is
true then return them in the standard coordinate system instead.
"""
function charge(
    u::StandaloneLightningEigenfunction{T},
    i::Integer,
    n::Integer,
    standard_coordinates = false,
) where {T<:Union{arb,fmpq}}
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
    k::Integer;
    boundary = nothing,
    notransform::Bool = false,
    n = div(length(coefficients(u)) - 1, 3) + 1, # Total number of charge points
) where {T <: Union{arb, arb_series}}
    if !notransform
        xy = coordinate_transformation(u, xy)
    end

    # Perform the change of coordinates corresponding to the charge
    # point
    xy = xy - charge(u, div(k - 1, 3) + 1, n)

    r = sqrt(xy[1]^2 + xy[2]^2)
    θ = atan(xy[2], xy[1])

    if mod1(k, 3) == 1
        return bessel_y(zero(λ), r*sqrt(λ))
    elseif mod1(k, 3) == 2
        return bessel_y(one(λ), r*sqrt(λ))*sin(θ)
    else
        return bessel_y(one(λ), r*sqrt(λ))*cos(θ)
    end
end

function (u::StandaloneLightningEigenfunction)(
    xy::AbstractVector{T},
    λ::arb,
    ks::UnitRange{Int};
    boundary = nothing,
    notransform::Bool = false,
    n = div(length(coefficients(u)) - 1, 3) + 1, # Total number of charge points
) where {T <: Union{arb, arb_series}}
    if !notransform
        xy = coordinate_transformation(u, xy)
    end

    res = similar(ks, T)
    i = 1
    charge_indices = (div(ks.start-1, 3)+1):(div(ks.stop-1, 3)+1)
    for charge_index in charge_indices
        # Perform the change of coordinates corresponding to the charge
        # point
        xy_local = xy - charge(u, charge_index, n)
        r, θ = polar_from_cartesian(xy_local)

        rsqrtλ = r*sqrt(λ)
        b = bessel_y(one(λ), rsqrtλ)
        s, c = sincos(θ)

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

    return res
end

# TODO: Figure out how to handle this
set_domain!(u::StandaloneLightningEigenfunction, ::AbstractDomain) = u

# TODO: Figure out how to handle this
recompute!(u::StandaloneLightningEigenfunction) = u
