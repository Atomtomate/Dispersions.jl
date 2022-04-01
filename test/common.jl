using FFTW
include("helper_functions.jl")

@testset "reduce_expand" begin
    for gr in grid_list
        for NN in 2:2:8
            kG = gen_kGrid(gr,4)
            arr1 = reduceKArr(kG, randn(rng, ComplexF64, gridshape(kG)))
            @test all(expandKArr(kG, arr1) .≈ expandKArr(kG, reduceKArr(kG, expandKArr(kG, arr1)))) # consistent?
            #TODO: test expand disp and direct disp calc
        end
    end
end

@testset "convolutions" begin
    for gr in grid_list
        kG = gen_kGrid(gr,4)
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
        arr1_sym = expandKArr(kG,reduceKArr(kG,arr1))
        arr2_sym = expandKArr(kG,reduceKArr(kG,arr2))
        conv_naive_res_2 = Dispersions.conv_sample_post(kG, naive_conv(arr1_sym, arr2_sym) ./ Nk(kG))
        conv_res = conv(kG, reduceKArr(kG, arr1_sym), reduceKArr(kG, arr2_sym))

        @testset "$gr" begin
            @test all(conv_theo_res .≈ conv_naive_fft_def_res)      # naive convolution and conv. theorem match
            @test all(conv_naive_res .≈ conv_theo_res_2)
            @test all(expandKArr(kG,conv_res) .≈ conv_naive_res_2)
            @test all(reduceKArr(kG, expandKArr(kG, conv_res)) .≈ conv_res)     # is symmetry preserved?
        end
    end
end
