@testset "KSum" begin
    for grid in grid_list
        for NN = 2:2:8
            gr = gen_kGrid(grid, NN)
            @test kintegrate(gr, ones(size(gr.ϵkGrid)))[1] ≈ 1.0
            @test_throws ArgumentError kintegrate(gr, ones(7, 18), 2)
            @test all(kintegrate(gr, ones(length(gr.ϵkGrid)..., 4), 1) .≈ ones(4))
        end
    end
end
