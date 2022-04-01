using Base.Iterators
#TODO: test kvices/ifft

@testset "gen_kGrid" begin
    gl = map(x -> gen_kGrid(x, NN), grid_list)
    for i = 1:length(gl)
        for NN in [4, 6, 8]
            kG = gl[i]
            Di = grid_list_D[i]
            @test kG.Nk == NN^Di
            @test kG.Ns == NN
            @test all(size(kG.fft_cache) .== repeat([NN], Di))
        end
    end
end
