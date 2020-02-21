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

export AbstractDomain,
    area, boundary_parameterization, boundary_points, interior_points,

    SphericalTriangle,
    angles, vertex, center, greatcircleplane, greatcircle,

    LShape,

    AbstractEigenfunction,
    coefficients, set_eigenfunction!,

    SphericalVertexEigenfunction, mu,

    SphericalInteriorEigenfunction,

    SphericalCombinedEigenfunction, basis_function,

    LShapeEigenfunction,

    sigma, sigma_coefficients,

    cartesian, spherical,

    mps, iteratemps,

    certifyindices

include("arb.jl")
include("arb_series.jl")
include("utilities.jl")

# Include types for domains and eigenfunctions
include("domain/domain_types.jl")
include("eigenfunction/eigenfunction_types.jl")

# Include methods for domains
include("domain/AbstractDomain.jl")
include("domain/SphericalTriangle.jl")
include("domain/LShape.jl")

# Include methods for eigenfunctions
include("eigenfunction/AbstractEigenfunction.jl")
include("eigenfunction/AbstractSphericalEigenfunction.jl")
include("eigenfunction/SphericalVertexEigenfunction.jl")
include("eigenfunction/SphericalInteriorEigenfunction.jl")
include("eigenfunction/SphericalCombinedEigenfunction.jl")
include("eigenfunction/LShapeEigenfunction.jl")

include("sigma.jl")
include("enclose.jl")
include("mps.jl")
include("index.jl")

end # module
