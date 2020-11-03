@recipe function f(u::StandaloneVertexEigenfunction{T}) where {T}
    seriestype := [:path :path :scatter]
    linewidth --> [2 3 3]
    markersize --> 5
    label --> ["" "u"]
    seriescolor --> [:blue :blue :red]
    linestyle --> [:dot :solid :solid]

    v = u.vertex
    orientation = ifelse(T == arb, u.orientation, u.parent(u.orientation)*u.parent(π))
    θ = ifelse(T == arb, u.θ, u.parent(u.θ)*u.parent(π))

    vertex = (Float64[u.vertex[1]], Float64[u.vertex[2]])
    edges = begin
        s1, c1 = sincos(orientation)
        v1 = v + 0.1*SVector(c1, s1)
        s2, c2 = sincos(orientation + θ)
        v2 = v + 0.1*SVector(c2, s2)
        (Float64[v1[1], u.vertex[1], v2[1]], Float64[v1[2], u.vertex[2], v2[2]])
    end
    arc = begin
        res = [
            v + 0.1*SVector(cos(orientation + t*θ), sin(orientation + t*θ))
            for t in range(0, 1, length = 20)
        ]
        (Float64.(getindex.(res, 1)), Float64.(getindex.(res, 2)))
    end

    x = hcat(
        arc[1],
        [edges[1]; fill(NaN, length(arc[1]) - length(edges[1]))],
        [vertex[1]; fill(NaN, length(arc[1]) - length(vertex[1]))],
    )
    y = hcat(
        arc[2],
        [edges[2]; fill(NaN, length(arc[2]) - length(edges[2]))],
        [vertex[2]; fill(NaN, length(arc[2]) - length(vertex[2]))],
    )

    return (x, y)
end
