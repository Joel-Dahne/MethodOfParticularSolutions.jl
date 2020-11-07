function eigenfunction_plotdata(u::StandaloneVertexEigenfunction{T}) where {T}
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

    return vertex, edges, arc
end

function eigenfunction_plotdata(u::StandaloneInteriorEigenfunction)
    vertex = (Float64[u.vertex[1]], Float64[u.vertex[2]])
    arc = begin
        res = [
            u.vertex + 0.1*SVector(cos(θ), sin(θ))
            for θ in range(0, 2π, length = 50)
        ]
        (Float64.(getindex.(res, 1)), Float64.(getindex.(res, 2)))
    end

    return vertex, arc
end

@recipe function f(u::StandaloneVertexEigenfunction)
    seriestype := [:path :path :scatter]
    linewidth --> [2 3 3]
    markersize --> 5
    label --> ["" "u"]
    seriescolor --> [:blue :blue :red]
    linestyle --> [:dot :solid :solid]

    vertex, edges, arc = eigenfunction_plotdata(u)

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

@recipe function f(u::StandaloneInteriorEigenfunction)
    seriestype := [:path :scatter]
    linewidth --> [2 3]
    markersize --> 5
    label --> ["" "u"]
    seriescolor --> [:blue :red]
    linestyle --> [:dot :solid]
    vertex, arc = eigenfunction_plotdata(u)

    x = hcat(
        arc[1],
        [vertex[1]; fill(NaN, length(arc[1]) - length(vertex[1]))],
    )
    y = hcat(
        arc[2],
        [vertex[2]; fill(NaN, length(arc[2]) - length(vertex[2]))],
    )

    return (x, y)
end

@userplot EigenfunctionPlot

@recipe function f(h::EigenfunctionPlot)
    if length(h.args) != 1 || !(typeof(h.args[1]) <: CombinedEigenfunction)
        error("EigenfunctionPlots should be given a CombinedEigenfunction. Got: $(typeof(h.args))")
    end
    u = h.args[1]

    # For some reason this makes it work... I don't understand plot
    # recipes...
    @series begin
        label := ""
        u.us[1]
    end

    for v in u.us
        @series begin
            label := ""
            v
        end
    end

    return nothing
end

@userplot EigenfunctionHeatmap

@recipe function f(h::EigenfunctionHeatmap)
    if !(length(h.args) == 3  || length(h.args) == 5) ||
        !(typeof(h.args[1]) <: AbstractDomain) ||
        !(typeof(h.args[2]) <: AbstractEigenfunction) ||
        !(typeof(h.args[3]) <: Union{Real,arb})
        error("EigenfunctionHeatmap should be given an eigenfunction, a domain and x, y. Got: $(typeof(h.args))")
    end
    if length(h.args) == 5
        domain, u, λ, xs, ys = h.args
    else
        domain, u, λ = h.args
        xs = 50
        ys = 50
    end

    vs = vertices(domain)
    if typeof(xs) <: Integer
        xs = range(extrema(Float64.(collect(getindex.(vs, 1))))..., length = xs)
    end
    if typeof(ys) <: Integer
        ys = range(extrema(Float64.(collect(getindex.(vs, 2))))..., length = ys)
    end

    pts = SVector.(domain.parent.(xs'), domain.parent.(ys));
    res = Float64.(abs.(u.(pts, domain.parent(λ))))
    @series begin
        seriestype := :heatmap
        Float64.(xs), Float64.(ys), res
    end

    @series begin
        label := ""
        domain, 100, 0
    end

    for v in u.us
        @series begin
            label := ""
            v
        end
    end

    return nothing
end