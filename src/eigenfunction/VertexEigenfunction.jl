function VertexEigenfunction(domain::Triangle, vertex::Integer; stride::Integer = 1)
    return VertexEigenfunction(
        domain,
        vertex,
        StandaloneVertexEigenfunction(domain, vertex, stride = stride),
    )
end

function Base.show(io::IO, u::VertexEigenfunction)
    println(io, "Vertex eigenfunction from vertex $(u.vertex)")
    if !haskey(io, :compact) || !io[:compact]
        println(io, "domain: $(u.domain)")
        print(io, "number of set coefficients: $(length(coefficients(u)))")
    end
end

coefficients(u::VertexEigenfunction) = u.u.coefficients
set_eigenfunction!(u::VertexEigenfunction, coefficients::Vector) =
    set_eigenfunction!(u.u, coefficients)

function active_boundaries(domain::Triangle, u::VertexEigenfunction)
    if domain === u.domain
        return u.vertex:u.vertex
    else
        return boundaries(domain)
    end
end

coordinate_transformation(u::VertexEigenfunction, xy::AbstractVector) =
    coordinate_transformation(u.u, xy)

function (u::VertexEigenfunction)(xy::AbstractVector{T},
                                  λ::arb,
                                  k::Integer;
                                  boundary = nothing,
                                  notransform::Bool = false,
                                  ) where {T <: Union{arb, arb_series}}
    return u.u(xy, λ, k, boundary = boundary, notransform = notransform)
end
