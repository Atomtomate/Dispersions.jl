using Base.Iterators

@testset "2D" begin
    sc_2d_2 = gen_cP_kGrid(2,2,1.0)
    sc_2d_16 = gen_cP_kGrid(16,2,1.0)
    @test Nk(sc_2d_2) == 2^2
    @test all(isapprox.(flatten(collect(gridPoints(sc_2d_2))), flatten([(0,0) (0,π); (π,0) (π,π)])))
    r2 = reduceKGrid(sc_2d_2)
    r16 = reduceKGrid(sc_2d_16)
    @test all(expandKGrid(r2, r2.ϵkGrid)[:] .≈ sc_2d_2.ϵkGrid)
    @test maximum(abs.(expandKGrid(r16, r16.ϵkGrid)[:] .- sc_2d_16.ϵkGrid)) < 1/10^12
end


@testset "3D" begin
    sc_3d_2 = gen_cP_kGrid(2, 3, 1.0)
    sc_3d_16 = gen_cP_kGrid(16, 3, 1.0)
    indTest = reshape([(1, 1, 1) (2, 1, 1) (1, 2, 1) (2, 2, 1) (1, 1, 2) (2, 1, 2) (1, 2, 2) (2, 2, 2)], (2,2,2))
    gridTest = reshape([(0, 0, 0) (π, 0, 0) (0, π, 0) (π, π, 0) (0, 0, π) (π, 0, π) (0, π, π) (π, π, π)], (2,2,2))
    @test Nk(sc_3d_2) == 2^3
    @test all(isapprox.(flatten(collect(gridPoints(sc_3d_2))), flatten(gridTest)))

    r2 = reduceKGrid(sc_3d_2)
    r16 = reduceKGrid(sc_3d_16)
    @test all(expandKGrid(r2, r2.ϵkGrid)[:] .≈ sc_3d_2.ϵkGrid)
    @test maximum(abs.(expandKGrid(r16, r16.ϵkGrid)[:] .- sc_3d_16.ϵkGrid)) < 1/10^12
end
