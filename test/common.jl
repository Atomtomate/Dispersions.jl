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
        fft_rek_2 = fft(reverse(reshape(rek_2,gridshape(kG))))[:]
        ifft_post_arr = randn(gridshape(kG))
        ifft_post_res1 = Dispersions.ifft_post(kG, ifft_post_arr)
        ifft_post_res2 = deepcopy(ifft_post_arr)
        Dispersions.ifft_post!(kG, ifft_post_res2)
        @test sum(abs.(ifft_post_res1 .- ifft_post_res2)) < 1e-8
        t1 = zeros(ComplexF64,length(kG.ϵkGrid))
        t2 = zeros(ComplexF64,length(kG.ϵkGrid))
        t3 = naive_conv(kG, rek, rek_2)[:]
        b1 = deepcopy(rek)
        b2 = deepcopy(rek_2)
        r1 = conv(kG, rek, rek_2)
        @test all(rek .== b1)
        @test all(rek_2 .== b2)
        conv!(kG, t1, rek, rek_2)
        @test sum(abs.(r1 .- t1)) < 1e-8
        @test sum(abs.(r1 .- t3)) < 1e-8        # fft conv .- naive conv
        r2 = conv_fft1(kG, rek, fft_rek_2)
        conv_fft1!(kG, t2, rek, fft_rek_2)
        @test sum(abs.(r1 .- r2)) < 1e-8        # fft conv .- 1 pre computed fft conv
        @test sum(abs.(r2 .- t2)) < 1e-8        # fft conv .- 1 pre computed fft conv (inplace)
        r2 = conv_fft(kG, fft_rek, fft_rek_2)
        conv_fft!(kG, t2, fft_rek, fft_rek_2)
        @test sum(abs.(r1 .- r2)) < 1e-8        # fft conv .- 2 pre computed fft conv
        @test sum(abs.(r2 .- t2)) < 1e-8
    end
    end
end
