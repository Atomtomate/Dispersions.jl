using Test

using Dispersions

@testset "gen_kGrid" begin
    t1_r = [(0.0,), (π,)]
    t1_i, t1 = gen_kGrid(2, 1, min=0, max=π, include_min=true)
    t2_r = [(0.0, 0.0) (0.0, π); (π, 0.0) (π, π)]
    t2_i, t2 = gen_kGrid(2, 2, min=0, max=π, include_min=true)
    #println(collect.(zip.(t2_r, (t2))))
    @test all(collect(t1_i) .== [(1,), (2,)])
    @test all([all(collect(t1)[i] .≈ t1_r[i]) for i in 1:length(t1)])
    @test all(collect(t2_i) .== [(1,1) (1,2); (2,1) (2,2)])
    @test all(all.(map(
                  x -> map(
                           y -> y[1] ≈ y[2],
                           x),
                  collect.(zip.(t2_r, t2))
                 )
             ))
end

