@testset "KSum" begin
    for NN in 1:8
        gr2 = gen_cP_kGrid(NN,2,1.3)
        gr3 = gen_cP_kGrid(NN,3,1.3)
        gr2_r = reduceKGrid(gr2)
        gr3_r = reduceKGrid(gr3)
        @test kintegrate(gr2, ones(size(gr2.ϵkGrid))) ≈ 1
        @test kintegrate(gr3, ones(size(gr3.ϵkGrid))) ≈ 1
        @test_throws ArgumentError kintegrate(gr2, ones(7,8))
        @test_throws ArgumentError kintegrate(gr2_r, ones(7,8))
        @test all(kintegrate(gr2_r, ones(length(gr2_r.ϵkGrid)..., 4), dim=1) .≈ ones(4))
        @test all(kintegrate(gr3_r, ones(length(gr3_r.ϵkGrid)..., 4), dim=1) .≈ ones(4))
    end
end
