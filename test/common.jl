using FFTW
# TODO error test for interface functions
kG_2 = gen_kGrid("2Dsc-1.4",16)
kG_3 = gen_kGrid("3Dsc-1.4",16)
hexG = gen_kGrid("p6m-1.4",16)
@testset "convolutions" begin
    for kG in map(x-> gen_kGrid(x,16), ["2Dsc-1.4", "3Dsc-1.4", "p6m-1.4"])
        fek = convert.(Complex{Float64},expandKArr(kG, kG.ϵkGrid))
        rek = convert.(Complex{Float64},deepcopy(kG.ϵkGrid))
        t1 = zeros(Complex{Float64},length(kG.ϵkGrid))
        t2 = zeros(Complex{Float64},length(kG.ϵkGrid))
        rek_2 = deepcopy(rek)
        fft_rek = fft(fek)
        r1 = conv(kG, rek, rek_2)
        conv!(kG, t1, rek, rek_2)
        @test all(r1 .≈ t1)
        r2 = conv_fft1(kG, rek, fft_rek)
        conv_fft1!(kG, t2, rek, fft_rek)
        @test all(r1 .≈ r2)
        @test all(r2 .≈ t2)
        r2 = conv_fft(kG, fft_rek, fft_rek)
        conv_fft!(kG, t2, fft_rek, fft_rek)
        @test all(r1 .≈ r2)
        @test all(r2 .≈ t2)
    end
end
