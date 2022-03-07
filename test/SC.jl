using Base.Iterators


@testset "2D" begin
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
    @test maximum(rr) < 1e-10
end


@testset "3D" begin
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
    @test maximum(rr) < 1e-10
end

@testset "reduce_expand" begin
    for NN in 1:16
        gl_gr, gl_w  = Dispersions.gausslegendre(NN);
        gr2 = Dispersions.gen_kGrid("2Dsc-1.3",NN, full=true)
        gr3 = Dispersions.gen_kGrid("3Dsc-1.3",NN, full=true)
        gr2_r = Dispersions.reduceKGrid(gr2)
        gr3_r = Dispersions.reduceKGrid(gr3)
        ek2 = reshape(gr2.ϵkGrid, (NN,NN))
        ek3 = reshape(gr3.ϵkGrid, (NN,NN,NN))
        res_r2 = zeros(eltype(gr2.ϵkGrid), size(gr2_r.ϵkGrid)...) 
        res_r3 = zeros(eltype(gr3.ϵkGrid), size(gr3_r.ϵkGrid)...) 
        res_f2 = zeros(eltype(ek2), size(ek2)...) 
        res_f3 = zeros(eltype(ek3), size(ek3)...) 
        @test all(reduceKArr(gr2_r, ek2) .≈ gr2_r.ϵkGrid)
        @test all(reduceKArr(gr3_r, ek3) .≈ gr3_r.ϵkGrid)
        reduceKArr!(gr2_r, res_r2, ek2)
        reduceKArr!(gr3_r, res_r3, ek3)
        @test all(res_r2 .≈ gr2_r.ϵkGrid)
        @test all(res_r3 .≈ gr3_r.ϵkGrid)
        gr2_cut = cut_mirror(ek2)
        gr3_cut = cut_mirror(ek3)
        @test all(map(x-> x in gr2.ϵkGrid, gr2_cut)) # No data lost
        @test all(map(x-> x in gr3.ϵkGrid, gr3_cut))
        al = ceil(Int,NN/2)
        gr2_pre_exp = zeros(size(ek2))
        gr2_pre_exp[al:end,al:end] = gr2_cut
        gr3_pre_exp = zeros(size(ek3))
        gr3_pre_exp[al:end,al:end,al:end] = gr3_cut
        @test all(abs.(expandKArr(gr2_r, gr2_r.ϵkGrid) .- ek2) .< 1.0/10^10)
        @test all(abs.(expandKArr(gr3_r, gr3_r.ϵkGrid) .- ek3) .< 1.0/10^10)
        expandKArr!(gr2_r, convert.(Complex{Float64},gr2_r.ϵkGrid))
        expandKArr!(gr3_r, convert.(Complex{Float64},gr3_r.ϵkGrid))
        @test all(abs.(gr2_r.expand_cache .- convert.(Complex{Float64},ek2)) .< 1.0/10^10)
        @test all(abs.(gr3_r.expand_cache .- convert.(Complex{Float64},ek3)) .< 1.0/10^10)
        #@test sum(gr2_r.kMult) == Nk(gr2_r)
        #@test sum(gr3_r.kMult) == Nk(gr3_r)
    end
end
