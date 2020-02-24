function triangle(i::Integer,
                  parent::ArbField = RealField(64);
                  withinterior = true)
    1 <= i <= 10 || throw(ArgumentError("give 1 ≤ i ≤ 10 not i = $i"))

    angles = [(3//4, 1//3, 1//2),
              (2//3, 1//3, 1//2),
              (2//3, 1//4, 1//2),
              (2//3, 1//3, 1//3),
              (3//4, 1//4, 1//3),
              (2//3, 1//4, 1//4),
              (2//3, 3//4, 3//4),
              (2//3, 2//3, 2//3),
              (1//2, 2//3, 3//4),
              (1//2, 2//3, 2//3),
              ]

    r = parent(1//8)
    λs = ball.(parent.(["12.400051652843377905",
                        "13.744355213213231835",
                        "20.571973537984730557",
                        "21.309407630190445260",
                        "24.456913796299111694",
                        "49.109945263284609920",
                        "4.261735",
                        "5.159146",
                        "6.241748",
                        "6.777108",
                        ]), r)

    domain = SphericalTriangle(fmpq.(angles[i]), parent)

    if i in 1:6
        stride = i in [4, 6] ? 2 : 1
        u = SphericalVertexEigenfunction(domain, 1, stride = stride)
    elseif i == 7
        if withinterior
            us = [SphericalVertexEigenfunction(domain, 1, stride = 2),
                  SphericalVertexEigenfunction(domain, 2, stride = 1),
                  SphericalVertexEigenfunction(domain, 3, stride = 1),
                  SphericalInteriorEigenfunction(domain, stride = 2)]
            u = SphericalCombinedEigenfunction(domain, us, [1, 1, 1, 4])
        else
            us = [SphericalVertexEigenfunction(domain, 1, stride = 2),
                  SphericalVertexEigenfunction(domain, 2, stride = 1),
                  SphericalVertexEigenfunction(domain, 3, stride = 1)]
            u = SphericalCombinedEigenfunction(domain, us, [1, 1, 1])
        end
    elseif i == 8
        if withinterior
            us = [[SphericalVertexEigenfunction(domain, i, stride = 2) for i in 1:3];
                  SphericalInteriorEigenfunction(domain, stride = 1)] #6
            u = SphericalCombinedEigenfunction(domain, us, [1, 1, 1, 4])
        else
            us = [[SphericalVertexEigenfunction(domain, i, stride = 2) for i in 1:3]]
            u = SphericalCombinedEigenfunction(domain, us, [1, 1, 1])
        end
    elseif i == 9
        if withinterior
            us = [SphericalVertexEigenfunction(domain, 2),
                  SphericalVertexEigenfunction(domain, 3),
                  SphericalInteriorEigenfunction(domain)]
            u = SphericalCombinedEigenfunction(domain, us, [1, 1, 4])
        else
            us = [SphericalVertexEigenfunction(domain, 2),
                  SphericalVertexEigenfunction(domain, 3)]
            u = SphericalCombinedEigenfunction(domain, us, [1, 1])
        end
    elseif i == 10
        if withinterior
            us = [SphericalVertexEigenfunction(domain, 2),
                  SphericalVertexEigenfunction(domain, 3),
                  SphericalInteriorEigenfunction(domain, stride = 2)]
            u = SphericalCombinedEigenfunction(domain, us, [1, 1, 4])
        else
            us = [SphericalVertexEigenfunction(domain, 2),
                  SphericalVertexEigenfunction(domain, 3)]
            u = SphericalCombinedEigenfunction(domain, us, [1, 1])
        end
    end

    (domain, u, λs[i])
end
