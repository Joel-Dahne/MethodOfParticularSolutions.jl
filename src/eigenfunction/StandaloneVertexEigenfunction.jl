function StandaloneVertexEigenfunction(
    vertex::SVector{2,arb},
    orientation::T,
    θ::T,
    stride::Integer = 1,
    parent::ArbField = parent(vertex[1]),
) where {T<:Union{arb,fmpq}}
    return StandaloneVertexEigenfunction(vertex, orientation, θ, stride, arb[], parent)
end

function StandaloneVertexEigenfunction(
    domain::Triangle{T},
    i::Integer;
    stride::Integer = 1,
    outside = false,
) where {T<:Union{arb,fmpq}}
    i ∈ boundaries(domain) ||
        throw(ArgumentError("attempt to get vertex $i from a $(typeof(domain))"))

    θ = domain.angles[i]

    if i == 1
        orientation = zero(θ)
    elseif i == 2
        orientation = ifelse(T == arb, domain.parent(π), 1) - θ
    elseif i == 3
        orientation = 2ifelse(T == arb, domain.parent(π), 1) - domain.angles[2] - θ
    end

    if outside
        orientation += θ
        θ = ifelse(T == arb, 2domain.parent(π), 2) - θ
    end

    return StandaloneVertexEigenfunction(
        vertex(domain, i),
        orientation,
        θ,
        stride,
        arb[],
        domain.parent,
    )
end

function StandaloneVertexEigenfunction(
    domain::TransformedDomain,
    i::Integer;
    stride::Integer = 1,
    outside = false,
)
    u = StandaloneVertexEigenfunction(
        domain.original,
        i,
        stride = stride,
        outside = outside,
    )
    # FIXME: This only works if u.θ::arb or if u.θ and
    # domain.orientation both are fmpq.
    u = StandaloneVertexEigenfunction(
        domain.map(u.vertex),
        u.orientation + domain.rotation,
        u.θ,
        stride,
        u.coefficients,
        domain.parent,
    )

    return u
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
nu(u::StandaloneVertexEigenfunction{fmpq}, k::Integer = 1) = u.parent(k * inv(u.θ))
nu(u::StandaloneVertexEigenfunction{arb}, k::Integer = 1) = u.parent(k * u.parent(π) / u.θ)

"""
    coordinate_transformation(u::StandaloneVertexEigenfunction, xy::AbstractVector)

Takes a 2-element vector `xy` representing a point in the plane in
Cartesian coordinates and makes a (affine) change of coordinates so
that `u.vertex` is put at the origin and rotated `u.orientation`
clockwise.
"""
function coordinate_transformation(
    u::StandaloneVertexEigenfunction{T},
    xy::AbstractVector,
) where {T<:Union{arb,fmpq}}
    if T == fmpq
        s, c = sincospi(-u.orientation, u.parent)
    elseif T == arb
        s, c = sincos(-u.orientation)
    end
    M = SMatrix{2,2}(c, s, -s, c)
    return M * (xy .- u.vertex)
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
        θ += 2parent(θ)(π)
    end

    rsqrtλ = r * sqrt(λ)
    res = similar(ks, T)
    for i in eachindex(ks)
        k = 1 + (ks[i] - 1) * u.stride
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
