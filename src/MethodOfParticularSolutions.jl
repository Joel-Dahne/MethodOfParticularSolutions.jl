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

using TimerOutputs
using Printf
using RecipesBase

export AbstractDomain,
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
include("special-functions.jl")

# Include types for domains and eigenfunctions
include("domain/types.jl")
include("domain/planar/types.jl")
include("domain/spherical/types.jl")
include("eigenfunction/eigenfunction_types.jl")

# Include methods for domains
include("domain/AbstractDomain.jl")
# Planar domains
include("domain/planar/AbstractPlanarDomain.jl")
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
include("eigenfunction/AbstractPlanarEigenfunction.jl")
include("eigenfunction/AbstractSphericalEigenfunction.jl")
include("eigenfunction/SphericalVertexEigenfunction.jl")
include("eigenfunction/SphericalInteriorEigenfunction.jl")
include("eigenfunction/SphericalCombinedEigenfunction.jl")
include("eigenfunction/KrewerasEigenfunction.jl")
include("eigenfunction/LShapeEigenfunction.jl")
include("eigenfunction/StandaloneVertexEigenfunction.jl")
include("eigenfunction/StandaloneInteriorEigenfunction.jl")
include("eigenfunction/StandaloneLightningEigenfunction.jl")
include("eigenfunction/LinkedEigenfunction.jl")
include("eigenfunction/CombinedEigenfunction.jl")
include("eigenfunction/plotting.jl")

include("sigma.jl")
include("enclose.jl")
include("mps.jl")
include("index.jl")
include("trace.jl")

include("examples/triangles.jl")
include("examples/denominators.jl")
include("examples/domains.jl")

end # module
