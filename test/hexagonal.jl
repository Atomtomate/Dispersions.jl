@testset "constructor" begin
    r2 = gen_kGrid("p6m-1.3",2)
    r3 = gen_kGrid("p6m-1.3",4)
    r16 = gen_kGrid("p6m-1.4",16)
    @test_throws ArgumentError expandKArr(r16, [1,2,3,4])
    #TODO: test gridpoints
    @test Nk(r2) == 2^2
    @test all(gridshape(r2) .== (2,2))
    @test all(dispersion(r2) .≈ r2.ϵkGrid)
end

@testset "kmult" begin
    kG = gen_kGrid("p6m-1.3",4)
    @test sum(kG.kMult) == Nk(kG)
end

@testset "reduce_expand" begin
    for NN in 2:5
        gr2 = Dispersions.gen_kGrid("p6m-1.3",NN, full=true)
        gr2_r = Dispersions.reduceKGrid(gr2)
        ek2 = reshape(gr2.ϵkGrid, (NN,NN))
        @test all(reduceKArr(gr2_r, ek2) .≈ gr2_r.ϵkGrid)
        gr2_cut = cut_mirror(ek2)
        @test all(map(x-> x in gr2.ϵkGrid, gr2_cut)) # No data lost
        @test all(abs.(expandKArr(gr2_r, gr2_r.ϵkGrid) .- ek2) .< 1.0/10^10)
    end
end
