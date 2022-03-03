using FFTW
include("helper_functions.jl")

@testset "general function" begin
    for NN in [4,6]
    gl = map(x-> gen_kGrid(x,NN), grid_list)
    for i in 1:length(gl)
        kG = gl[i]
        @test Nk(kG) == kG.Ns^grid_list_D[i]
        @test gridshape(kG) == Tuple(repeat([kG.Ns], grid_list_D[i]))
    end
    end
end

# TODO error test for interface functions
# TODO: this needs loads of tests for different lattices
@testset "convolutions" begin
    for NN in [3,4,5]
    for kG in map(x-> gen_kGrid(x,NN), grid_list)
        rek = convert.(ComplexF64,deepcopy(kG.ϵkGrid))
        rek_2 = randn(ComplexF64, size(rek))
        fft_rek = fft(reshape(rek,gridshape(kG)))[:]
        ifft_post_arr = randn(gridshape(kG))
        ifft_post_res1 = Dispersions.ifft_post(kG, ifft_post_arr)
        ifft_post_res2 = similar(ifft_post_arr)
        Dispersions.ifft_post!(kG, ifft_post_res2, ifft_post_arr)
        @test all(ifft_post_res1 .≈ ifft_post_res2)
        t1 = zeros(ComplexF64,length(kG.ϵkGrid))
        t3 = naive_conv(kG, rek, rek_2)[:]
        r1 = conv(kG, rek, rek_2)
        t2 = conv_old(kG, rek, rek_2)
        conv!(kG, t1, rek, rek_2)
        @test all(r1 .≈ t1)
        @test all(r1 .≈ t2)
        @test abs(sum(r1 .- t3)) < 10e-8
        r2 = conv_fft1(kG, rek, fft_rek)
        conv_fft1!(kG, t2, rek, fft_rek)
        @test abs(sum(r1 .- r2)) < 10e-8
        @test all(r2 .≈ t2)
        r2 = conv_fft(kG, fft_rek, fft_rek)
        conv_fft!(kG, t2, fft_rek, fft_rek)
        @test abs(sum(r1 .- r2)) < 10e-8
        @test all(r2 .≈ t2)
    end
    end
end
