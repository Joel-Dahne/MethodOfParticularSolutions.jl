function example_domain_H(parent = RealField(precision(BigFloat)))
    angles = fmpq.([3//2, 1//2, 1//2, 1//2, 1//2, 3//2, 3//2, 1//2, 1//2, 1//2, 1//2, 3//2])
    vertices = [parent.(v) for v in [
        (0, 0),
        (0, 1),
        (-1, 1),
        (-1, -2),
        (0, -2),
        (0, -1),
        (1, -1),
        (1, -2),
        (2, -2),
        (2, 1),
        (1, 1),
        (1, 0),
    ]]

    domain = Polygon(angles, vertices, parent)

    u1 = StandaloneVertexEigenfunction(vertex(domain, 1), fmpq(1//2), fmpq(3//2))
    u2 = StandaloneVertexEigenfunction(vertex(domain, 6), fmpq(0), fmpq(3//2))
    u3 = StandaloneVertexEigenfunction(vertex(domain, 7), fmpq(-1//2), fmpq(3//2))
    u4 = StandaloneVertexEigenfunction(vertex(domain, 12), fmpq(1), fmpq(3//2))


    u = CombinedEigenfunction(
        domain,
        [u1, u2, u3, u4],
        us_to_boundary = [
            setdiff(boundaries(domain), (12, 1)),
            setdiff(boundaries(domain), (5, 6)),
            setdiff(boundaries(domain), (6, 7)),
            setdiff(boundaries(domain), (11, 12))
        ],
    )

    return domain, u
end

function example_domain_ngon(n, parent = RealField(precision(BigFloat)))
    θ = fmpq((n - 2)//n)
    angles = fill(θ, n)
    vertices = [(cos(θ), sin(θ)) for θ in (2parent(π)/n).*(0:n-1)]
    domain = Polygon(angles, vertices, parent)

    us = [
        StandaloneVertexEigenfunction(vertex(domain, i), i*(1 - θ) + θ*1//2, θ)
        for i in boundaries(domain)
    ]

    u = CombinedEigenfunction(
        domain,
        us,
        us_to_boundary = [
            setdiff(boundaries(domain), (i, mod1(i - 1, length(boundaries(domain)))))
            for i in boundaries(domain)
        ]
    )

    return domain, u
end
