using Base.Iterators
#TODO: test kvices/ifft

@testset "gen_kGrid" begin
    for NN in [3,4,5]
    gl = map(x-> gen_kGrid(x,NN), grid_list)
    for i in 1:length(gl)
        kG = gl[i]
        Di = grid_list_D[i]
        @test kG.Nk == NN^Di
        @test kG.Ns == NN
        @test all(size(kG.fft_cache) .== repeat([NN],Di))
    end
    end
end
