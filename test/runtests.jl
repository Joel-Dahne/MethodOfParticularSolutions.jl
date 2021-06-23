using MethodOfParticularSolutions, Nemo, Printf, Test

@testset "MethodOfParticularSolutions" begin
    include("sphericaltriangle.jl")
    include("lshape.jl")
    include("triangles.jl")
    include("example-domains.jl")
end
