@testset "spherical triangles" begin
    Ns = 8:8:16

    RR = RealField(1024)
    results = RR.(["[12.400051652843377905286053412896636720736 +/- 1.64e-40]",
                   "[13.7443552132132318354011215921380207828066502596318748941363320689579830 +/- 3.27e-71]",
                   "[20.5719735379847305566258422 +/- 5.38e-26]",
                   "[21.3094076301904452589534814412305177783368425771467166131131424182062385470402339419123020595676115778838298367063775989397269169412254133009366735802749167865869428407055350499081173154930 +/- 4.29e-188]",
                   "[24.45691379629911169448043814477268289960795913156636922934413916 +/- 7.30e-63]",
                   "[49.10994526328460991967034315150826835369842561533395606847954650063727524833998848617655899444520661743928451538721837069883497076 +/- 4.90e-129]",
                   "[4.26173475529398708575 +/- 7.62e-21]",
                   "[5.159145642466541711221675 +/- 3.96e-25]",
                   "[6.24174833072633424 +/- 7.59e-18]",
                   "[6.7771080545983009573857 +/- 5.76e-23]",
                   ])

    goalradius = [2e-5, 4e-9, 5e-4, 3e-20, 2e-7, 2e-12, 1.2e-1, 2e-1, 2e-1, 9e-4]

    for i in 1:10
        domain, u, interval = MethodOfParticularSolutions.triangle(i)
        println("$i: Computing eigenvalue for the $domain")

        optim_prec_final = ifelse(i != 4, 60, 100)
        λs = iteratemps(domain, u, interval, Ns,
                        optim_prec_final = optim_prec_final,
                        optim_prec_linear = true,
                        show_trace = true)
        rad = Float64(radius(λs[end]))

        @test overlaps(results[i], λs[end])
        @test rad < goalradius[i]

        if rad < goalradius[i]
            @printf "radius = %e < %e\n" rad goalradius[i]
        else
            @printf "radius = %e >= %e\n" rad goalradius[i]
        end
        println("")
    end
end
