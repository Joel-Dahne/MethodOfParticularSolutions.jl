module MPS

using Nemo
using Random
using CoordinateTransformations
using LinearAlgebra
using GenericSVD
using Optim

export SphericalTriangle,
    angles, vertex, greatcircleplane, greatcircle,
    boundary_points, interior_points,

    set_eigenfunction!,

    SphericalVertexEigenfunction, mu,

    sigma, sigma_coefficients,

    cartesian, spherical,

    mps

include("arb.jl")
include("utilities.jl")

include("domain.jl")
include("eigenfunction.jl")

include("SphericalTriangle.jl")

include("SphericalVertexEigenfunction.jl")

include("sigma.jl")
include("mps.jl")

end # module
