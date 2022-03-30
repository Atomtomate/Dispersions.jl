using FFTW
include("helper_functions.jl")
# TODO error test for interface functions
# TODO: this needs loads of tests for different lattices
@testset "convolutions" begin
    for gr in grid_list
        kG = gen_kGrid(gr,4)
        arr1 = randn(rng, gridshape(kG))
        arr2 = randn(rng, gridshape(kG))
        arr1_fft = fft(arr1)
        arr2_fft = fft(arr2)
        conv_theo_res = ifft(arr1_fft .* arr2_fft)
        conv_naive_res = naive_conv_fft_def(arr1, arr2)
        @testset "$gr" begin
            @test all(conv_theo_res .≈ conv_naive_res)
        end
    end
end
