@testset "KSum" begin
    for grid in grid_list
        for NN in 2:8
            gr = gen_kGrid(grid,NN)
            @test kintegrate(gr, ones(size(gr.ϵkGrid)))[1] ≈ 1.0
            @test_throws ArgumentError kintegrate(gr, ones(7,8))
            @test all(kintegrate(gr, ones(length(gr.ϵkGrid)..., 4), dim=1) .≈ ones(4))
        end
    end
end
