using Base.Iterators


@testset "2D" begin
    # @test_throws ArgumentError gen_kGrid("2Dsc-1.3",3)
    r2 = gen_kGrid("2Dsc-1.3--1.4-1.5",2)
    @test r2.t ≈ 1.3
    @test r2.tp ≈ -1.4
    @test r2.tpp ≈ 1.5
    r16 = gen_kGrid("2Dsc-1.4--1.5-1.6",16)
    @test Nk(r2) == 2^2
    @test all(dispersion(r2) .≈ r2.ϵkGrid)
    @test_throws ArgumentError expandKArr(r16, [1,2,3,4])
    @test all(gridshape(r2) .== (2,2))
    @test isapprox(kintegrate(r16, r16.ϵkGrid), 0.0, atol=1e-10)
end
