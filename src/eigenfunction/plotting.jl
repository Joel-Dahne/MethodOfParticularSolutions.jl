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

function eigenfunction_plotdata(u::StandaloneInteriorEigenfunction{T}) where {T}
    vertex = (Float64[u.vertex[1]], Float64[u.vertex[2]])
    edge = begin
        orientation = ifelse(T == arb, u.orientation, u.parent(u.orientation)*u.parent(π))
        s, c = sincos(orientation)
        v = u.vertex + 0.1*SVector(c, s)
        (Float64[u.vertex[1], v[1]], Float64[u.vertex[2], v[2]])
    end
    arc = begin
        res = [
            u.vertex + 0.1*SVector(cos(θ), sin(θ))
            for θ in range(0, 2π, length = 50)
        ]
        (Float64.(getindex.(res, 1)), Float64.(getindex.(res, 2)))
    end

    return vertex, edge, arc
end

function eigenfunction_plotdata(u::StandaloneLightningEigenfunction{T}) where {T}
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

    charges = begin
        n = div(length(coefficients(u)) - 1, 3) + 1
        cs = [charge(u, i, n, true) for i in 1:n]
        (Float64.(getindex.(cs, 1)), Float64.(getindex.(cs, 2)))
    end

    return vertex, edges, arc, charges
end

@recipe function f(u::StandaloneVertexEigenfunction)
    seriestype := [:path :path :scatter]
    linewidth := [2 3 3]
    markersize := 5
    label --> ["" "u" ""]
    seriescolor := [:blue :blue :red]
    linestyle := [:dot :solid :solid]

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
    seriestype := [:path :path :scatter]
    linewidth := [2 2 3]
    markersize := 5
    label --> ["" "" "u"]
    seriescolor := [:blue :blue :red]
    linestyle := [:dot :solid :solid]
    vertex, edge, arc = eigenfunction_plotdata(u)

    x = hcat(
        arc[1],
        [edge[1]; fill(NaN, length(arc[1]) - length(edge[1]))],
        [vertex[1]; fill(NaN, length(arc[1]) - length(vertex[1]))],
    )
    y = hcat(
        arc[2],
        [edge[2]; fill(NaN, length(arc[2]) - length(edge[2]))],
        [vertex[2]; fill(NaN, length(arc[2]) - length(vertex[2]))],
    )
    return (x, y)
end

@recipe function f(u::StandaloneLightningEigenfunction)
    seriestype --> [:path :path :scatter :scatter]
    linewidth := [2 3 3 2]
    markersize := [5 5 5 2]
    label --> ["" "" "u" ""]
    seriescolor := [:blue :blue :red :green]
    linestyle := [:dot :solid :solid :solid]

    vertex, edges, arc, charges = eigenfunction_plotdata(u)

    x = hcat(
        arc[1],
        [edges[1]; fill(NaN, length(arc[1]) - length(edges[1]))],
        [vertex[1]; fill(NaN, length(arc[1]) - length(vertex[1]))],
        [charges[1]; fill(NaN, length(arc[1]) - length(charges[1]))],
    )
    y = hcat(
        arc[2],
        [edges[2]; fill(NaN, length(arc[2]) - length(edges[2]))],
        [vertex[2]; fill(NaN, length(arc[2]) - length(vertex[2]))],
        [charges[2]; fill(NaN, length(arc[2]) - length(charges[2]))],
    )

    return (x, y)
end

@userplot EigenfunctionPlot

@recipe function f(h::EigenfunctionPlot)
    if length(h.args) != 1 || !(typeof(h.args[1]) <: CombinedEigenfunction)
        error("EigenfunctionPlots should be given a CombinedEigenfunction. Got: $(typeof(h.args))")
    end
    u = h.args[1]

    for v in u.us
        @series begin
            label := ""
            v
        end
    end

    return nothing
end

@userplot EigenfunctionHeatmap

@recipe function f(
    h::EigenfunctionHeatmap;
    absolute_value = false,
    twosided = !absolute_value,
    include_exterior = false,
)
    if !(length(h.args) == 3  || length(h.args) == 5) ||
        !(typeof(h.args[1]) <: AbstractDomain) ||
        !(typeof(h.args[2]) <: AbstractEigenfunction) ||
        !(typeof(h.args[3]) <: Union{Real,arb})
        throw(ArgumentError(
            "EigenfunctionHeatmap should be given a domain, an eigenfunction, an eigenvalue " *
            "and x, y. Got: $(typeof(h.args))"
        ))
    end
    if length(h.args) == 5
        domain, u, λ, xs, ys = h.args[1:5]
    elseif length(h.args) == 3
        domain, u, λ = h.args[1:3]
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
    res = similar(pts, Float64)
    let λ = domain.parent(λ)
        @Threads.threads for i in eachindex(pts)
            if include_exterior || pts[i] ∈ domain
                res[i] = ifelse(absolute_value, abs, identity)(u(pts[i], λ))
            else
                res[i] = 0
            end
        end
    end

    xlims := extrema(xs)
    ylims := extrema(ys)
    label := ""
    aspect_ratio --> :equal

    if twosided
        # Make the color bar symmetric around 0
        m = maximum(abs, res)
        clims --> (-m, m)
    end

    @series begin
        seriestype := :heatmap
        if twosided
            seriescolor --> :balance
        else
            seriescolor --> :matter
        end
        xs, ys, res
    end

    @series begin
        domain, 100, 0
    end

    return nothing
end
