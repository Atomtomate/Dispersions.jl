using Base.Iterators

@testset "dispersion" begin
    r8 = gen_kGrid("fcc-1.2",2)
    r64 = gen_kGrid("fcc-1.1",4)
    indTest = reduceKArr(r8, reshape([(1, 1, 1) (2, 1, 1) (1, 2, 1) (2, 2, 1) (1, 1, 2) (2, 1, 2) (1, 2, 2) (2, 2, 2)], (2,2,2)))
    gridTest = reduceKArr(r8, reshape([(0, 0, 0) (-π, π, π) (π, -π, π) (0, 0, 2π) (π, π, -π) (0, 2π, 0) (2π, 0, 0) (π, π, π)], (2,2,2)))
    @test Nk(r8) == 2^3
    @test Nk(r64) == 4^3
    @test all(dispersion(r8) .≈ r8.ϵkGrid)
    @test all(isapprox.(flatten(gridPoints(r8)), flatten(gridTest)))
    @test_throws ArgumentError expandKArr(r64, [1,2,3,4])
    @test all(gridshape(r8) .== (2,2,2))
    @test isapprox(kintegrate(r64, r64.ϵkGrid), 0.0, atol=num_eps)
    @test isapprox(kintegrate(r64, r64.ϵkGrid .* r64.ϵkGrid), 6 * r64.t^2, atol=num_eps)
    #rr = abs.(conv(r64, convert.(ComplexF64,r64.ϵkGrid), convert.(ComplexF64,r64.ϵkGrid)) .- ( - r64.t .* r64.ϵkGrid))
    #@test maximum(rr) < 1e-10
end
