@testset "KSum" begin
    for NN in 2:2
        gr2 = gen_kGrid("2Dsc-1.3",NN)
        gr3 = gen_kGrid("3Dsc-1.3",NN)
        gr2_r = reduceKGrid(gr2)
        gr3_r = reduceKGrid(gr3)
        @test kintegrate(gr2_r, ones(size(gr2_r.ϵkGrid)))[1] ≈ 1.0
        @test kintegrate(gr3_r, ones(size(gr3_r.ϵkGrid)))[1] ≈ 1.0
        @test_throws ArgumentError kintegrate(gr2_r, ones(7,8))
        @test all(kintegrate(gr2_r, ones(length(gr2_r.ϵkGrid)..., 4), dim=1) .≈ ones(4))
        @test all(kintegrate(gr3_r, ones(length(gr3_r.ϵkGrid)..., 4), dim=1) .≈ ones(4))
    end
end
