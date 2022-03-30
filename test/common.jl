using FFTW
include("helper_functions.jl")
# TODO error test for interface functions
# TODO: this needs loads of tests for different lattices
@testset "convolutions" begin
    for kG in map(x-> gen_kGrid(x,4),grid_list)
        
        fek = convert.(ComplexF64,expandKArr(kG, kG.ϵkGrid))
        rek = convert.(ComplexF64,deepcopy(kG.ϵkGrid))
        rek_2 = deepcopy(rek)
        t1 = zeros(ComplexF64,length(kG.ϵkGrid))
        t2 = zeros(ComplexF64,length(kG.ϵkGrid))
        test_data = randn(ComplexF64, length(kG.kMult))
        test2_data = randn(ComplexF64, length(kG.kMult))
        test2_data_expanded = expandKArr(kG, test2_data)
        test_data_expanded = expandKArr(kG, test_data)

        conv_naive_data_expanded = naive_conv(kG, rek, rek_2)
        conv_naive_data = reduceKArr(kG, conv_naive_data_expanded)

        fft_rek = fft(fek)
        r1 = conv(kG, rek, rek_2)
        conv!(kG, t1, rek, rek_2)
        @test all(r1 .≈ t1)
        @test all(abs.(conv_naive_data .- r1) .< num_eps)
        r2 = conv_fft1(kG, rek, fft_rek)
        conv_fft1!(kG, t2, rek, fft_rek)
        @test all(r1 .≈ r2)
        @test all(r2 .≈ t2)
        r2 = conv_fft(kG, fft_rek, fft_rek)
        conv_fft!(kG, t2, fft_rek, fft_rek)
        @test all(r1 .≈ r2)
        @test all(r2 .≈ t2)

        
        #@test all(conv_naive_data .≈ r2)
        #@test all(conv_naive_data .≈ t2)
        #=
        rek = convert.(ComplexF64,expandKArr(kG, kG.ϵkGrid)[:])
        rek_2 = randn(ComplexF64, size(rek))
        ones_2 = ones(ComplexF64, size(rek))
        rek_red = convert.(ComplexF64,deepcopy(kG.ϵkGrid)[:])
        rek_2_red = randn(ComplexF64, size(rek))
        ones_2_red = ones(ComplexF64, size(rek))
        println(size(rek))
        #fft_rek = fft(reshape(rek,gridshape(kG)))[:]
        #fft_rek_2 = fft(reverse(reshape(rek_2,gridshape(kG))))[:]
        #ifft_post_arr = randn(gridshape(kG))
        #ifft_post_res1 = Dispersions.ifft_post(kG, ifft_post_arr)
        #ifft_post_res2 = deepcopy(ifft_post_arr)
        #Dispersions.ifft_post!(kG, ifft_post_res2)
        #@test sum(abs.(ifft_post_res1 .- ifft_post_res2)) < 1e-8
        t1 = zeros(ComplexF64,length(kG.ϵkGrid))
        t2 = naive_conv(kG, rek, rek_2)[:]
        t3 = conv(kG, rek, ones_2)
        t4 = conv(kG, ones_2, rek)
        b1 = deepcopy(rek)
        b2 = deepcopy(rek_2)
        r1 = conv(kG, rek, rek_2)
        @test all(rek .== b1)
        @test all(rek_2 .== b2)
        conv!(kG, t1, rek, rek_2)
        @test sum(abs.(r1 .- t1)) < 1e-8
        @test sum(abs.(r1 .- t2)) < 1e-8        # fft conv .- naive conv
        @test all(abs.(sum(rek)/Nk(kG) .- t3) .< 1e-8 )
        @test all(abs.(sum(rek)/Nk(kG) .- t4) .< 1e-8 )
        r2 = conv_fft1(kG, rek, fft_rek_2)
        conv_fft1!(kG, t2, rek, fft_rek_2)
        @test sum(abs.(r1 .- r2)) < 1e-8        # fft conv .- 1 pre computed fft conv
        @test sum(abs.(r2 .- t2)) < 1e-8        # fft conv .- 1 pre computed fft conv (inplace)
        r2 = conv_fft(kG, fft_rek, fft_rek_2)
        conv_fft!(kG, t2, fft_rek, fft_rek_2)
        @test sum(abs.(r1 .- r2)) < 1e-8        # fft conv .- 2 pre computed fft conv
        @test sum(abs.(r2 .- t2)) < 1e-8
        =#
    end
end
