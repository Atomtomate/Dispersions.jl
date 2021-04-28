@testset "KSum" begin
    for NN in 2:3
        gr2 = gen_cP_kGrid(NN,2,1.3)
        gr3 = gen_cP_kGrid(NN,3,1.3)
        gr2_r = reduceKGrid(gr2)
        gr3_r = reduceKGrid(gr3)
        @test kintegrate(gr2, ones(size(gr2.ϵkGrid))) ≈ 1
        @test kintegrate(gr3, ones(size(gr3.ϵkGrid))) ≈ 1
    end
end
