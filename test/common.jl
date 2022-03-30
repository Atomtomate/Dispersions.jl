using FFTW
include("helper_functions.jl")
# TODO error test for interface functions
# TODO: this needs loads of tests for different lattices
@testset "convolutions" begin
    for gr in grid_list
        kG = gen_kGrid(gr,3)
        arr1 = randn(rng, ComplexF64, gridshape(kG))
        arr2 = randn(rng, ComplexF64, gridshape(kG))
        arr1_fft = fft(arr1)
        arr1_r_fft = fft(reverse(arr1))
        arr2_r_fft = fft(reverse(arr2))
        arr2_fft = fft(arr2)
        conv_theo_res = ifft(arr1_fft .* arr2_fft)
        conv_naive_fft_def_res = naive_conv_fft_def(arr1, arr2)     # ∑_k arr1[k] * arr2[q-k]
        conv_naive_res = naive_conv(arr1, arr2)                     # ∑_k arr1[k] * arr2[q+k]
        conv_theo_res_2 = reverse(ifft(arr1_fft .* arr2_r_fft))
        conv_theo_res_3 = reverse(ifft(arr1_r_fft .* arr2_fft))
        conv_res = conv(kG, reduceKArr(kG, arr1), reduceKArr(kG, arr2))
        conv_naive_res_2 = naive_conv(expandKArr(kG,reduceKArr(kG,arr1)), expandKArr(kG,reduceKArr(kG,arr2)))

        @testset "$gr" begin
            @test all(conv_theo_res .≈ conv_naive_fft_def_res)      # naive convolution and conv. theorem match
            @test all(conv_naive_res .≈ conv_theo_res_2)
            #@test all(conv_theo_res_3 .≈ conv_theo_res_2)
            @test all(expandKArr(kG,conv_res) .≈ conv_naive_res_2 ./ Nk(kG))
        end
    end
end
