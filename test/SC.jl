using Base.Iterators

@testset "2D" begin
    sc_2d_2 = gen_cP_kGrid(2,2,1.3)
    sc_2d_3 = gen_cP_kGrid(4,2,1.3)
    sc_2d_16 = gen_cP_kGrid(16,2,1.4)
    @test Nk(sc_2d_2) == 2^2
    @test all(isapprox.(flatten(collect(gridPoints(sc_2d_2))), flatten([(0,0) (0,π); (π,0) (π,π)])))
    r2 = reduceKGrid(sc_2d_2)
    r3 = reduceKGrid(sc_2d_2)
    r16 = reduceKGrid(sc_2d_16)
    #@test all(expandKGrid(r2, r2.ϵkGrid)[:] .≈ sc_2d_2.ϵkGrid)
    #@test all(expandKGrid(r3, r3.ϵkGrid)[:] .≈ sc_2d_3.ϵkGrid)
    #@test maximum(abs.(expandKGrid(r16, r16.ϵkGrid)[:] .- sc_2d_16.ϵkGrid)) < 1/10^12
end


@testset "3D" begin
    sc_3d_2 = gen_cP_kGrid(2, 3, 1.2)
    sc_3d_16 = gen_cP_kGrid(4, 3, 1.1)
    indTest = reshape([(1, 1, 1) (2, 1, 1) (1, 2, 1) (2, 2, 1) (1, 1, 2) (2, 1, 2) (1, 2, 2) (2, 2, 2)], (2,2,2))
    gridTest = reshape([(0, 0, 0) (π, 0, 0) (0, π, 0) (π, π, 0) (0, 0, π) (π, 0, π) (0, π, π) (π, π, π)], (2,2,2))
    @test Nk(sc_3d_2) == 2^3
    @test all(isapprox.(flatten(collect(gridPoints(sc_3d_2))), flatten(gridTest)))

    r2 = reduceKGrid(sc_3d_2)
    r16 = reduceKGrid(sc_3d_16)
    #@test all(expandKGrid(r2, r2.ϵkGrid)[:] .≈ sc_3d_2.ϵkGrid)
    #@test maximum(abs.(expandKGrid(r16, r16.ϵkGrid)[:] .- sc_3d_16.ϵkGrid)) < 1/10^12
end

@testset "internal" begin
    # mirror
    for NN in 2:16
        gr2 = gen_cP_kGrid(NN,2,1.3)
        gr3 = gen_cP_kGrid(NN,3,1.3)
        gr2_r = reduceKGrid(gr2)
        gr3_r = reduceKGrid(gr3)
        ek2 = reshape(gr2.ϵkGrid, (NN,NN))
        ek3 = reshape(gr3.ϵkGrid, (NN,NN,NN))
        gr2_cut = Dispersions.cut_mirror(ek2)
        gr3_cut = Dispersions.cut_mirror(ek3)
        @test all(map(x-> x in gr2.ϵkGrid, gr2_cut)) # No data lost
        @test all(map(x-> x in gr3.ϵkGrid, gr3_cut))
        @test all(abs.(Dispersions.expand_mirror(gr2_r.kInd, gr2_cut) .- ek2) .< 1.0/10^10)
    end
end
