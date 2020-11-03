function StandaloneVertexEigenfunction(
    vertex::SVector{2,arb},
    orientation::T,
    θ::T,
    stride::Integer = 1,
    parent::ArbField = parent(vertex[1]),
) where {T <: Union{arb,fmpq}}
    return StandaloneVertexEigenfunction(
        vertex,
        orientation,
        θ,
        stride,
        arb[],
        parent,
    )
end

function StandaloneVertexEigenfunction(
    domain::Triangle{T},
    i::Integer;
    stride::Integer = 1,
) where {T <: Union{arb,fmpq}}
    i ∈ boundaries(domain) || throw(ArgumentError("attempt to get vertex $i from a $(typeof(domain))"))

    θ = domain.angles[i]
    if i == 1
        orientation = zero(θ)
    elseif i == 2
        orientation = ifelse(T == arb, domain.parent(π), 1) - θ
    elseif i == 3
        orientation = 2ifelse(T == arb, domain.parent(π), 1) -
            domain.angles[2] - θ
    end

    return StandaloneVertexEigenfunction(
        vertex(domain, i),
        orientation,
        domain.angles[i],
        stride,
        arb[],
        domain.parent,
    )
end

function Base.show(io::IO, u::StandaloneVertexEigenfunction)
    println(io, "Standalone Vertex eigenfunction")
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
nu(u::StandaloneVertexEigenfunction{fmpq}, k::Integer = 1) = u.parent(-k*inv(u.θ))
nu(u::StandaloneVertexEigenfunction{arb}, k::Integer = 1) = u.parent(-k*u.parent(π)/u.θ)

"""
        coordinate_transformation(u::StandaloneVertexEigenfunction, xy::AbstractVector)

    Takes a 2-element vector `xy` representing a point in the plane in
    Cartesian coordinates and makes a (affine) change of coordinates so
    that `u.vertex` is put at the origin and rotated `u.orientation`
    clockwise.
    """
function coordinate_transformation(u::StandaloneVertexEigenfunction{fmpq}, xy::AbstractVector)
    s, c = sincospi(-u.orientation, u.parent)
    M = SMatrix{2, 2}(c, s, -s, c)
    return M*(xy .- u.vertex)
end

function coordinate_transformation(u::StandaloneVertexEigenfunction{arb}, xy::AbstractVector)
    s, c = sincos(-u.orientation)
    M = SMatrix{2, 2}(c, s, -s, c)
    return M*(xy .- u.vertex)
end

"""
        coordinate_transformation(u::StandaloneVertexEigenfunction, r, θ)

    Takes a point `r, θ` in polar coordinates (affine) change of
    coordinates so the vertex `u` originates from is put at the origin
    with the right edge on the x-axis.
    """
function coordinate_transformation(u::StandaloneVertexEigenfunction, r, θ)
    return polar_from_cartesian(coordinate_transformation(u, cartesian_from_polar(r, θ)))
end

function (u::StandaloneVertexEigenfunction)(xy::AbstractVector{T},
                                            λ::arb,
                                            k::Integer;
                                            boundary = nothing,
                                            notransform::Bool = false,
                                            ) where {T <: Union{arb, arb_series}}
    if !notransform
        xy = coordinate_transformation(u, xy)
    end

    k = 1 + (k - 1)*u.stride

    ν = nu(u, k)
    r, θ = polar_from_cartesian(xy)
    return bessel_j(ν, sqrt(λ)*r)*sin(ν*θ)
end

function (u::StandaloneVertexEigenfunction)(r::T,
                                            θ::T,
                                            λ::arb,
                                            k::Integer;
                                            boundary = nothing,
                                            notransform::Bool = false,
                                            ) where {T <: Union{arb, arb_series}}
    if !notransform
        r, θ = coordinate_transformation(u, r, θ)
    end

    k = 1 + (k - 1)*u.stride

    ν = nu(u, k)
    return bessel_j(ν, sqrt(λ)*r)*sin(ν*θ)
end

# TODO: Figure out how to handle this
set_domain!(u::StandaloneVertexEigenfunction, _::AbstractDomain) = u

# TODO: Figure out how to handle this
recompute!(u::StandaloneVertexEigenfunction) = u

# Have a function that computes u from the domain?
