using MethodOfParticularSolutions, Nemo, Printf, Test

@testset "MethodOfParticularSolutions" begin
    @testset "end-to-end-tests" begin
        include("end-to-end-tests/sphericaltriangle.jl")
        include("end-to-end-tests/lshape.jl")
        include("end-to-end-tests/triangles.jl")
        include("end-to-end-tests/example-domains.jl")
    end
end
