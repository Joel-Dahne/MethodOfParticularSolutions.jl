function VertexEigenfunction(domain::Triangle, vertex::Integer; stride::Integer = 1)
    vertex ∈ boundaries(domain) || throw(ArgumentError("attempt to get vertex $vertex from a $(typeof(domain))"))
    return VertexEigenfunction(domain, vertex, stride, arb[])
end

function Base.show(io::IO, u::VertexEigenfunction)
    println(io, "Vertex eigenfunction from vertex $(u.vertex)")
    if !haskey(io, :compact) || !io[:compact]
        println(io, "domain: $(u.domain)")
        print(io, "number of set coefficients: $(length(u.coefficients))")
    end
end

function active_boundaries(u::VertexEigenfunction)
    u.vertex .== boundaries(u.domain)
end

"""
    nu(u::VertexEigenfunction, k::Integer)

Return `k*ν₀` as an `arb`, the parameter used for the Bessel function.
"""
function nu(u::VertexEigenfunction{fmpq}, k::Integer = 1)
    u.domain.parent(-k*inv(angledivπ(u.domain, u.vertex)))
end

"""
    coordinate_transformation(u::VertexEigenfunction, xy::AbstractVector)

Takes a 2-element vector `xy` representing a point in the plane in
Cartesian coordinates and makes a (affine) change of coordinates so
the vertex `u` originates from is put at the origin with the right
edge on the x-axis.
"""
function coordinate_transformation(u::VertexEigenfunction, xy::AbstractVector)
    if u.vertex == 1
        return xy
    end
    if u.vertex == 2
        θ = angle(u.domain, 2) - u.domain.parent(π)
    elseif u.vertex == 3
        θ = angle(u.domain, 2) + angle(u.domain, 3) - 2u.domain.parent(π)
    end
    s, c = sincos(θ)
    M = SMatrix{2, 2}(c, s, -s, c)
    return M*(xy .- vertex(u.domain, u.vertex))
end

"""
    coordinate_transformation(u::VertexEigenfunction, r, θ)

Takes a point `r, θ` in polar coordinates (affine) change of
coordinates so the vertex `u` originates from is put at the origin
with the right edge on the x-axis.
"""
function coordinate_transformation(u::VertexEigenfunction, r, θ)
    return polar_from_cartesian(coordinate_transformation(u, cartesian_from_polar(r, θ)))
end

function (u::VertexEigenfunction)(xy::AbstractVector{T},
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

function (u::VertexEigenfunction)(r::T,
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
