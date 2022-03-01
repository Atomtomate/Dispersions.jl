using FFTW
include("helper_functions.jl")
# TODO error test for interface functions
# TODO: this needs loads of tests for different lattices
@testset "convolutions" begin
    for NN in [4,6]
    for kG in map(x-> gen_kGrid(x,NN), ["2Dsc-1.4", "3Dsc-1.4", "p6m-1.4", "fcc-1.4"])
        fek = convert.(ComplexF64,expandKArr(kG, kG.ϵkGrid))
        rek = convert.(ComplexF64,deepcopy(kG.ϵkGrid))
        t1 = zeros(ComplexF64,length(kG.ϵkGrid))
        t2 = zeros(ComplexF64,length(kG.ϵkGrid))
        test_data = randn(ComplexF64, length(kG.kMult))
        test2_data = randn(ComplexF64, length(kG.kMult))
        test2_data_expanded = expandKArr(kG, test2_data)
        test_data_expanded = expandKArr(kG, test_data)

        conv_naive_data_expanded = naive_conv(kG, test_data, test2_data)
        conv_naive_data = reduceKArr(kG, conv_naive_data_expanded)

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
end
