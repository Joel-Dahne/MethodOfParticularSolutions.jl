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

function example_domain_ngon(
    n,
    parent = RealField(precision(BigFloat));
    lightning = false,
    linked = false,
)
    θ = fmpq((n - 2)//n)
    angles = fill(θ, n)
    vertices = [(cos(θ), sin(θ)) for θ in (2parent(π)/n).*(0:n-1)]
    domain = Polygon(angles, vertices, parent)

    if !lightning
        us = [
            StandaloneVertexEigenfunction(vertex(domain, i), i*(1 - θ) + θ*1//2, θ)
            for i in boundaries(domain)
        ]
        us_to_boundary = [
            setdiff(boundaries(domain), (i, mod1(i - 1, length(boundaries(domain)))))
            for i in eachindex(us)
        ]
    else
        if linked
            us = [
                LinkedEigenfunction(
                    [
                        StandaloneLightningEigenfunction(vertex(domain, i), i*(1 - θ) + θ*1//2, θ)
                        for i in boundaries(domain)
                    ],
                    [1, 1, 1, 1],
                ),
                StandaloneInteriorEigenfunction(domain),
            ]
        else
            us = [[
                StandaloneLightningEigenfunction(vertex(domain, i), i*(1 - θ) + θ*1//2, θ)
                for i in boundaries(domain)
            ]; StandaloneInteriorEigenfunction(domain)]
        end
        us_to_boundary = fill(boundaries(domain), length(us))
    end

    u = CombinedEigenfunction(domain, us, us_to_boundary = us_to_boundary)

    return domain, u
end

function example_domain_triangle_in_triangle(parent = RealField(precision(BigFloat)))
    domain1 = Triangle(fmpq(1//3), fmpq(1//3), parent)
    domain2 = TransformedDomain(
        Triangle(fmpq(1//3), fmpq(1//3), parent),
        fmpq(0),
        parent(0.25),
        SVector(parent(0.375), parent(0.375))
    )

    domain = IntersectedDomain(domain1, domain2)

    us = vcat(
        [StandaloneLightningEigenfunction(domain1, i) for i in 1:3],
        [StandaloneLightningEigenfunction(domain2, i, outside = true, l = parent(0.15)) for i in 1:3],
        [StandaloneInteriorEigenfunction([parent(0.5), parent(0.375/2)])]
    )

    u = CombinedEigenfunction(
        domain,
        us,
        [1, 1, 1, 1, 1, 1, 1],
        us_to_boundary = [
            1:6,
            1:6,
            1:6,
            1:6,
            1:6,
            1:6,
            1:6,
        ]
    )

    return domain, u
end

function example_domain_ngon_in_ngon(n1, n2, parent = RealField(precision(BigFloat)))
    θ1 = fmpq((n1 - 2)//n1)
    angles1 = fill(θ1, n1)
    vertices1 = [(cos(θ), sin(θ)) for θ in (2parent(π)/n1).*(0:n1-1)]
    domain1 = Polygon(angles1, vertices1, parent)

    θ2 = fmpq((n2 - 2)//n2)
    angles2 = fill(θ2, n2)
    vertices2 = [(cos(θ), sin(θ)) for θ in (2parent(π)/n2).*(0:n2-1)]
    domain2 = TransformedDomain(
        Polygon(angles2, vertices2, parent),
        fmpq(0),
        parent(0.5),
        SVector(parent(0), parent(0)),
    )

    domain = IntersectedDomain(domain1, domain2)

    us = vcat(
        [StandaloneLightningEigenfunction(domain1, i) for i in boundaries(domain1)],
        [
            StandaloneLightningEigenfunction(domain2, i, outside = true, l = parent(0.1))
            for i in boundaries(domain2)
        ],
        [StandaloneInteriorEigenfunction(domain)]
    )

    u = CombinedEigenfunction(
        domain,
        us,
        [1, 1, 1, 1, 1, 1, 1],
    )

    return domain, u
end

function example_domain_goal(parent = RealField(precision(BigFloat)))
    # The main domain is a hexagon
    n = 6
    θ = fmpq((n - 2)//n)
    angles = fill(θ, n)
    vertices = [(cospi(θ, parent), sinpi(θ, parent)) for θ in fmpq(2//n).*(0:n-1)]
    exterior = Polygon(angles, vertices, parent)

    # We remove six triangles
    domain_triangle = TransformedDomain(
        Triangle(fmpq(1//3), fmpq(1//3), parent),
        fmpq(0),
        parent(1//4),
        SVector(cospi(fmpq(1//6), parent), sinpi(fmpq(1//6), parent))/4,
    )

    interiors = [
        TransformedDomain(
            Triangle(fmpq(1//3), fmpq(1//3), parent),
            fmpq(i//3),
            parent(0.25),
            SVector(cospi(fmpq(i//3 + 1//6), parent), sinpi(fmpq(i//3 + 1//6), parent))/4
        )
        for i in 0:5
    ]

    domain = IntersectedDomain(exterior, interiors)

    # TODO: These should all have the same coefficients and be
    # symmetric (so we should skip the sin ones?)
    us1 = [
        StandaloneLightningEigenfunction(
            vertex(exterior, i),
            fmpq(mod(Rational(1 - θ//2 + (i - 1)*(1 - θ)), 2)),
            θ,
        )
        for i in boundaries(exterior)
    ]

    # TODO: Are there any symmetries here?
    us2 = [
        StandaloneLightningEigenfunction(d, 1, l = parent(0.15), outside = true)
        for d in interiors
    ]

    # TODO: Are there any symmetries here?
    us3 = [
        StandaloneLightningEigenfunction(d, i, l = parent(0.15), outside = true)
        for d in interiors, i in 2:3
    ][:]

    # TODO: Are there any symmetries here?
    us4 = [
        StandaloneInteriorEigenfunction(domain)
    ]

    u = CombinedEigenfunction(
        domain,
        vcat(us1, us2, us3, us4),
    )

    return domain, u
end
