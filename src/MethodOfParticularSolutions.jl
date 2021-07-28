module MethodOfParticularSolutions

using Nemo
using ArbTools
using Random
using StaticArrays
using CoordinateTransformations
using LinearAlgebra
using GenericSVD
using Optim
using OrderedCollections
import SpecialFunctions
import SpecialFunctions: besselj, besselj0, besselj1, bessely, bessely0, bessely1

using TimerOutputs
using Printf
using RecipesBase

export AbstractDomain,
    vertexindices,
    boundaries,
    area,
    boundary_parameterization,
    boundary_points,
    interior_points,
    SphericalTriangle,
    angles,
    vertex,
    center,
    greatcircleplane,
    greatcircle,
    LShape,
    Triangle,
    Polygon,
    TransformedDomain,
    IntersectedDomain,
    AbstractEigenfunction,
    coefficients,
    set_eigenfunction!,
    SphericalVertexEigenfunction,
    mu,
    SphericalInteriorEigenfunction,
    SphericalCombinedEigenfunction,
    basis_function,
    KrewerasEigenfunction,
    LShapeEigenfunction,
    StandaloneVertexEigenfunction,
    StandaloneInteriorEigenfunction,
    StandaloneLightningEigenfunction,
    LinkedEigenfunction,
    CombinedEigenfunction,
    sigma,
    sigma_coefficients,
    cartesian,
    spherical,
    mps!,
    iteratemps,
    certifyindices

include("arb.jl")
include("arb_series.jl")
include("utilities.jl")

# Include types for domains and eigenfunctions
include("domain/types.jl")
include("domain/planar/types.jl")
include("domain/spherical/types.jl")
include("eigenfunction/types.jl")
include("eigenfunction/planar/types.jl")
include("eigenfunction/spherical/types.jl")

# Include methods for domains
include("domain/AbstractDomain.jl")
# Planar domains
include("domain/planar/Triangle.jl")
include("domain/planar/Polygon.jl")
include("domain/planar/LShape.jl")
include("domain/planar/TransformedDomain.jl")
include("domain/planar/IntersectedDomain.jl")
include("domain/planar/plotting.jl")
# Spherical domains
include("domain/spherical/SphericalTriangle.jl")

# Include methods for eigenfunctions
include("eigenfunction/AbstractEigenfunction.jl")
# Planar eigenfunctions
include("eigenfunction/planar/AbstractPlanarEigenfunction.jl")
include("eigenfunction/planar/StandaloneVertexEigenfunction.jl")
include("eigenfunction/planar/StandaloneInteriorEigenfunction.jl")
include("eigenfunction/planar/StandaloneLightningEigenfunction.jl")
include("eigenfunction/planar/LShapeEigenfunction.jl")
include("eigenfunction/planar/LinkedEigenfunction.jl")
include("eigenfunction/planar/CombinedEigenfunction.jl")
include("eigenfunction/planar/plotting.jl")
# Spherical eigenfunctions
include("eigenfunction/spherical/AbstractSphericalEigenfunction.jl")
include("eigenfunction/spherical/SphericalVertexEigenfunction.jl")
include("eigenfunction/spherical/SphericalInteriorEigenfunction.jl")
include("eigenfunction/spherical/SphericalCombinedEigenfunction.jl")
include("eigenfunction/spherical/KrewerasEigenfunction.jl")

include("sigma.jl")
include("enclose.jl")
include("mps.jl")
include("index.jl")
include("trace.jl")

include("examples/triangles.jl")
include("examples/denominators.jl")
include("examples/domains.jl")

end # module
