@testset "regular triangles" begin
    N = 16
    prec = 64

    RR = RealField(prec)

    old_prec = precision(BigFloat)
    setprecision(BigFloat, prec)

    results = RR.(["12.40005165284337790528605341 +/- 4.96e-27",
                   "13.7443552132132318354011215921380207828066502596318748941363320689579830254389619211598920192317 +/- 7.13e-95",
                   "20.5719735379847305566258421533 +/- 6.17e-29",
                   "21.309407630190445258953481441230517778336842577146716613113142418206238547040233941912302059567611577883829836706377598939726916941225413300936673580274916786587 +/- 2.21e-160",
                   "24.4569137962991116944804381447726828996080 +/- 6.73e-41",
                   "49.109945263284609919670343151508268353698425615333956068479546500637275248339988486176558994445206617439284515387218370698834970763 +/- 3.91e-130"])

    for i in 1:6
        domain, interval = triangle(i, RR)
        u = SphericalVertexEigenfunction(domain, 1)

        λ, u = mps(domain, u, interval, N)

        @test overlaps(results[i], λ)
        @show results[i]
        @show λ
        @show Float64(radius(λ))
    end

    setprecision(BigFloat, old_prec)
end
