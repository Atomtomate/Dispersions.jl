@testset "KSum" begin
    for NN in [3,20,30,50,80]
    for grid in grid_list
            kG = gen_kGrid(grid,NN)
            #TODO: this hoilds for SC, find tests for all lattice types
            #TODO: for SC: @test isapprox(kintegrate(kG, kG.ϵkGrid .* kG.ϵkGrid), 6 * kG.t^2, atol=1e-10)
            @test isapprox(kintegrate(kG, kG.ϵkGrid), 0.0, atol=1e-10)
            @test kintegrate(kG, ones(size(kG.ϵkGrid)))[1] ≈ 1.0
            #@test_throws ArgumentError kintegrate(kG, ones(7,18), 2)
            @test all(kintegrate(kG, ones(length(kG.ϵkGrid)..., 4), 1) .≈ ones(4))
        end
    end
end
