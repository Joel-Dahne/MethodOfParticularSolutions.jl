module MPS

using Nemo
using Random
using CoordinateTransformations
using LinearAlgebra
using GenericSVD
using Optim

export SphericalTriangle,
    angles, vertex, center, greatcircleplane, greatcircle,
    boundary_points, interior_points,

    set_eigenfunction!,

    SphericalVertexEigenfunction, mu,

    SphericalInteriorEigenfunction,

    sigma, sigma_coefficients,

    cartesian, spherical,

    mps

include("arb.jl")
include("utilities.jl")

# Include types for domains and eigenfunctions
include("domain/domain_types.jl")
include("eigenfunction/eigenfunction_types.jl")

# Include methods for domains
include("domain/SphericalTriangle.jl")

# Include methods for eigenfunctions
include("eigenfunction/AbstractSphericalEigenfunction.jl")
include("eigenfunction/SphericalVertexEigenfunction.jl")
include("eigenfunction/SphericalInteriorEigenfunction.jl")

include("sigma.jl")
include("mps.jl")

end # module
