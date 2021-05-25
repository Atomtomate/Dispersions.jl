@testset "KSum" begin
    for grid in grid_list
        for NN in 2:2
            gr = gen_kGrid(grid,NN)
            gr_r = reduceKGrid(gr)
            @test kintegrate(gr_r, ones(size(gr_r.ϵkGrid)))[1] ≈ 1.0
            @test_throws ArgumentError kintegrate(gr_r, ones(7,8))
            @test all(kintegrate(gr_r, ones(length(gr_r.ϵkGrid)..., 4), dim=1) .≈ ones(4))
        end
    end
end
