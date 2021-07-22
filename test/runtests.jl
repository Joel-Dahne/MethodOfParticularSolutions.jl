using MethodOfParticularSolutions, Nemo, Printf, Test

@testset "MethodOfParticularSolutions" begin
    @testset "unit-tests" begin
        include("unit-tests/Triangle.jl")
        include("unit-tests/Polygon.jl")
        include("unit-tests/TransformedDomain.jl")
        include("unit-tests/IntersectedDomain.jl")
    end

    @testset "end-to-end-tests" begin
        include("end-to-end-tests/sphericaltriangle.jl")
        include("end-to-end-tests/lshape.jl")
        include("end-to-end-tests/triangles.jl")
        include("end-to-end-tests/example-domains.jl")
    end
end
