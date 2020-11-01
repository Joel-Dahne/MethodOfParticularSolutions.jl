function SphericalVertexEigenfunction(domain::SphericalTriangle,
                                      vertex::Int;
                                      stride::Int = 1)
    vertex >= 1 && vertex <= 3 || throw(ErrorException("vertex must be between 1 and 3 not $vertex"))
    SphericalVertexEigenfunction(domain, vertex, stride, arb[])
end

function Base.show(io::IO, u::SphericalVertexEigenfunction)
    println(io, "Spherical Vertex eigenfunction from vertex $(u.vertex)")
    if !haskey(io, :compact) || !io[:compact]
        println(io, "domain: $(u.domain)")
        print(io, "number of set coefficients: $(length(u.coefficients))")
    end
end

function active_boundaries(u::SphericalVertexEigenfunction)
    u.vertex .== (1, 2, 3)
end

"""
    mu(eigenfunction::SphericalVertexEigenfunction,
       k::Integer = 1)
> Return k*μ0 as an arb ball, the parameter used for the Legendre
  function.
"""
function mu(u::SphericalVertexEigenfunction{fmpq},
            k::Integer = 1)
    u.domain.parent(-k*inv(u.domain.angles[u.vertex]))
end

function mu(u::SphericalVertexEigenfunction{arb},
            k::Integer = 1)
    u.domain.parent(-k*inv(u.domain.angles[u.vertex]))*u.domain.parent(π)
end

function coordinate_transformation(u::SphericalVertexEigenfunction,
                                   xyz::AbstractVector{T}
                                   ) where {T <: Union{arb, arb_series}}
    # TODO: The performance could most likely be improved
    if u.vertex == 1
        return xyz
    elseif u.vertex == 2
        # Rotate by α along the y-axis so that the second vertex ends
        # up at the north pole, then rotate by β around the z-axis so
        # that the boundary ends up parallel to the y-axis.
        α = -acos(vertex(u.domain, 2)[3])
        β = angle(u.domain, 2) - u.domain.parent(π)
        L = LinearMap(RotZY(β, α))
        return L(xyz)
    elseif u.vertex == 3
        # Rotate by α along the z-axis so that the third vertex
        # ends up parallel to the y-axis, then rotate by β around
        # the y-axis so that it ends up at the north pole. Finally
        # rotate by γ so along the z-axis so that the boundary
        # ends up parallel to the y-axis.
        α = -angle(u.domain, 1)
        β = -acos(vertex(u.domain, 3)[3])
        γ = u.domain.parent(π)
        L = LinearMap(RotZYZ(γ, β, α))
        return L(xyz)
    else
        throw(ErrorException("vertex must be between 1 and 3 not $vertex"))
    end
end

function coordinate_transformation(u::SphericalVertexEigenfunction,
                                   θ::T,
                                   ϕ::T
                                   ) where {T <: Union{arb, arb_series}}
    if u.vertex == 1
        return θ, ϕ
    else
        return spherical(coordinate_transformation(u, cartesian(θ, ϕ)))
    end
end

function (u::SphericalVertexEigenfunction)(xyz::AbstractVector{T},
                                           λ::arb,
                                           k::Integer;
                                           boundary = nothing,
                                           notransform::Bool = false
                                           ) where {T <: Union{arb, arb_series}}
    if !isnothing(boundary) && boundary != u.vertex
        if T == arb
            return u.domain.parent(0)
        else
            return 0*xyz[1]
        end
    end
    k = 1 + (k - 1)*u.stride

    ν::arb = -0.5 + sqrt(0.25 + λ)
    μ::arb = mu(u, k)

    if !notransform
        xyz = coordinate_transformation(u, xyz)
    end

    ϕ = atan(xyz[2], xyz[1])
    legendre_p_safe(ν, μ, xyz[3])*sin(μ*ϕ)
end

function (u::SphericalVertexEigenfunction)(θ::T,
                                           ϕ::T,
                                           λ::arb,
                                           k::Integer;
                                           boundary = nothing,
                                           notransform::Bool = false
                                           ) where {T <: Union{arb, arb_series}}
    if !isnothing(boundary) && boundary != u.vertex
        if T == arb
            return u.domain.parent(0)
        else
            return 0*θ
        end
    end
    k = 1 + (k - 1)*u.stride

    ν::arb = -0.5 + sqrt(0.25 + λ)
    μ::arb = mu(u, k)

    if !notransform
        θ, ϕ = coordinate_transformation(u, θ, ϕ)
    end

    legendre_p_safe(ν, μ, cos(θ))*sin(μ*ϕ)
end

function norm(domain::SphericalTriangle,
              u::SphericalVertexEigenfunction,
              λ::arb)
    ϕ_integral = -domain.parent(π)/(2*mu(u, 1))

    CC = ComplexField(domain.parent.prec)
    # The integrals goes from zero to the lower bound for θ. However
    # the function has a branch cut at zero and Arb has problem
    # handling this. We therefore integrate a small distance away from
    # zero.
    a = CC(1e-1)

    # Create a new domain which has the vertex we are expanding from
    # on the north pole.
    domain = SphericalTriangle(domain.angles[[u.vertex,
                                              mod1(u.vertex + 1, 3),
                                              mod1(u.vertex + 2, 3)]],
                               domain.parent)
    b = CC(theta_bound(domain))

    θ_integral = domain.parent(0)
    for k in 1:min(4, length(u.coefficients))
        ν::arb = -0.5 + sqrt(0.25 + λ)
        μ::arb = mu(u, 1 + (k - 1)*u.stride)
        c2 = u.coefficients[k]^2

        f = θ -> begin
            θ = real(θ)
            if !isnonpositive(cos(θ) - 1)
                return CC(NaN)
            end

            CC(c2*sin(θ)legendre_p_safe(ν, μ, cos(θ))^2)
        end
        y = Nemo.integrate(CC, f, a, b)
        θ_integral += real(y)
    end

    sqrt(ϕ_integral*θ_integral)
end
