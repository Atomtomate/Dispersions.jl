using Base.Iterators

@testset "3D" begin
    r8 = gen_kGrid("fcc-1.2",2)
    r64 = gen_kGrid("fcc-1.1",4)
    indTest = reduceKArr(r8, reshape([(1, 1, 1) (2, 1, 1) (1, 2, 1) (2, 2, 1) (1, 1, 2) (2, 1, 2) (1, 2, 2) (2, 2, 2)], (2,2,2)))
    gridTest = reduceKArr(r8, reshape([(0, 0, 0) (π, 0, 0) (0, π, 0) (π, π, 0) (0, 0, π) (π, 0, π) (0, π, π) (π, π, π)], (2,2,2)))
    @test Nk(r8) == 2^3
    @test Nk(r64) == 4^3
    @test all(dispersion(r8) .≈ r8.ϵkGrid)
    @test all(isapprox.(flatten(gridPoints(r8)), flatten(gridTest)))
    @test_throws ArgumentError expandKArr(r64, [1,2,3,4])
    @test all(gridshape(r8) .== (2,2,2))
    @test isapprox(kintegrate(r64, r64.ϵkGrid), 0.0, atol=1e-10)
    @test isapprox(kintegrate(r64, r64.ϵkGrid .* r64.ϵkGrid), 6 * r64.t^2, atol=1e-10)
    rr = abs.(conv(r64, convert.(ComplexF64,r64.ϵkGrid), convert.(ComplexF64,r64.ϵkGrid)) .- ( - r64.t .* r64.ϵkGrid))
    @test maximum(rr) < 1e-10
end

@testset "reduce_expand" begin
    for NN in 1:16
        gr3 = Dispersions.gen_kGrid("3Dsc-1.3",NN, full=true)
        gr3_r = Dispersions.reduceKGrid(gr3)
        ek3 = reshape(gr3.ϵkGrid, (NN,NN,NN))
        res_r3 = zeros(eltype(gr3.ϵkGrid), size(gr3_r.ϵkGrid)...) 
        res_f3 = zeros(eltype(ek3), size(ek3)...) 
        @test all(reduceKArr(gr3_r, ek3) .≈ gr3_r.ϵkGrid)
        reduceKArr!(gr3_r, res_r3, ek3)
        @test all(res_r3 .≈ gr3_r.ϵkGrid)
        gr3_cut = cut_mirror(ek3)
        @test all(map(x-> x in gr3.ϵkGrid, gr3_cut))
        al = ceil(Int,NN/2)
        gr3_pre_exp = zeros(size(ek3))
        gr3_pre_exp[al:end,al:end,al:end] = gr3_cut
        @test all(abs.(expandKArr(gr3_r, gr3_r.ϵkGrid) .- ek3) .< 1.0/10^10)
        expandKArr!(gr3_r, convert.(Complex{Float64},gr3_r.ϵkGrid))
        @test all(abs.(gr3_r.expand_cache .- convert.(Complex{Float64},ek3)) .< 1.0/10^10)
    end
end
