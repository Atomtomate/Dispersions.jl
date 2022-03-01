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
#@testset "convolutions" begin
#    for NN in [4,6]
#    for kG in map(x-> gen_kGrid(x,NN), grid_list)
#        rek = convert.(ComplexF64,deepcopy(kG.ϵkGrid))
#        t1 = zeros(ComplexF64,length(kG.ϵkGrid))
#        t2 = zeros(ComplexF64,length(kG.ϵkGrid))
#        rek_2 = deepcopy(rek)
#        #TODO: naive conv
#        r1 = conv(kG, rek, rek_2)
#        conv!(kG, t1, rek, rek_2)
#        @test all(r1 .≈ t1)
#        r2 = conv_fft1(kG, rek, fft_rek)
#        conv_fft1!(kG, t2, rek, fft_rek)
#        @test all(r1 .≈ r2)
#        @test all(r2 .≈ t2)
#        r2 = conv_fft(kG, fft_rek, fft_rek)
#        conv_fft!(kG, t2, fft_rek, fft_rek)
#        @test all(r1 .≈ r2)
#        @test all(r2 .≈ t2)
#    end
#    end
#end
