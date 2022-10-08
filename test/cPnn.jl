using Base.Iterators


@testset "2D" begin
    @test_throws ArgumentError gen_kGrid("2Dsc-1.3",3)
    r2 = gen_kGrid("2Dsc-1.3",2)
    r16 = gen_kGrid("2Dsc-1.4",16)
    @test Nk(r2) == 2^2
    @test all(dispersion(r2) .≈ r2.ϵkGrid)
    @test all(isapprox.(flatten(gridPoints(r2)), flatten([(0,0) (π,0) (π,π)])))
    @test all(isapprox.(flatten(expandKArr(r2, gridPoints(r2))), flatten([(0,0) (π,0); (π,0) (π,π)])))
    @test_throws ArgumentError expandKArr(r16, [1,2,3,4])
    @test all(gridshape(r2) .== (2,2))
    @test isapprox(kintegrate(r16, r16.ϵkGrid), 0.0, atol=1e-10)
    @test isapprox(kintegrate(r16, r16.ϵkGrid .* r16.ϵkGrid), 4 * r16.t^2, atol=1e-10)
    rr = abs.(conv(r16, convert.(ComplexF64,r16.ϵkGrid), convert.(ComplexF64,r16.ϵkGrid)) .- ( - r16.t .* r16.ϵkGrid))
    #@test maximum(rr) < 1e-10
end


@testset "3D" begin
    @test_throws ArgumentError gen_kGrid("2Dsc-1.3",3)
    r2 = gen_kGrid("3Dsc-1.2",2)
    r16 = gen_kGrid("3Dsc-1.1",4)
    indTest = reduceKArr(r2, reshape([(1, 1, 1) (2, 1, 1) (1, 2, 1) (2, 2, 1) (1, 1, 2) (2, 1, 2) (1, 2, 2) (2, 2, 2)], (2,2,2)))
    gridTest = reduceKArr(r2, reshape([(0, 0, 0) (π, 0, 0) (0, π, 0) (π, π, 0) (0, 0, π) (π, 0, π) (0, π, π) (π, π, π)], (2,2,2)))
    @test Nk(r2) == 2^3
    @test Nk(r16) == 4^3
    @test all(dispersion(r2) .≈ r2.ϵkGrid)
    @test all(isapprox.(flatten(gridPoints(r2)), flatten(gridTest)))
    @test_throws ArgumentError expandKArr(r16, [1,2,3,4])
    @test all(gridshape(r2) .== (2,2,2))
    @test isapprox(kintegrate(r16, r16.ϵkGrid), 0.0, atol=1e-10)
    @test isapprox(kintegrate(r16, r16.ϵkGrid .* r16.ϵkGrid), 6 * r16.t^2, atol=1e-10)
    rr = abs.(conv(r16, convert.(ComplexF64,r16.ϵkGrid), convert.(ComplexF64,r16.ϵkGrid)) .- ( - r16.t .* r16.ϵkGrid))
    #@test maximum(rr) < 1e-10
end
