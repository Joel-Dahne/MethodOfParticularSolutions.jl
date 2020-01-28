module MPS

using Nemo
using ArbToolsNemo
using Random
using StaticArrays
using CoordinateTransformations
using LinearAlgebra
using GenericSVD
using Optim

using TimerOutputs

export SphericalTriangle,
    angles, vertex, center, greatcircleplane, greatcircle,
    boundary_parameterization, boundary_points, interior_points,

    coefficients, set_eigenfunction!,

    SphericalVertexEigenfunction, mu,

    SphericalInteriorEigenfunction,

    SphericalCombinedEigenfunction, basis_function,

    sigma, sigma_coefficients,

    cartesian, spherical,

    mps

include("arb.jl")
include("arb_series.jl")
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
include("eigenfunction/SphericalCombinedEigenfunction.jl")

include("sigma.jl")
include("enclose.jl")
include("mps.jl")

end # module
