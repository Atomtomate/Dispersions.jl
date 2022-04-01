@testset "constructor" begin
    r2 = gen_kGrid("p6m-1.3", 2)
    r3 = gen_kGrid("p6m-1.3", 4)
    r16 = gen_kGrid("p6m-1.4", 16)
    @test_throws ArgumentError expandKArr(r16, [1, 2, 3, 4])
    #TODO: test gridpoints
    @test Nk(r2) == 2^2
    @test all(gridshape(r2) .== (2, 2))
    @test all(dispersion(r2) .≈ r2.ϵkGrid)
end

@testset "kmult" begin
    kG = gen_kGrid("p6m-1.3", 4)
    @test sum(kG.kMult) == Nk(kG)
end
