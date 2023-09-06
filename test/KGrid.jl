using Base.Iterators
# TODO: test kvices/ifft
# TODO: test ϵ_k_plus_q for other grid types than simple cubic

@testset "gen_kGrid" begin
    @test false skip = true
    # gl = map(x -> gen_kGrid(x, NN), grid_list)
    # for i = 1:length(gl)
    #     for NN in [4, 6, 8]
    #         kG = gl[i]
    #         Di = grid_list_D[i]
    #         @test kG.Nk == NN^Di
    #         @test kG.Ns == NN
    #         @test all(size(kG.fft_cache) .== repeat([NN], Di))
    #     end
    # end
end

@testset "shifted grid" begin
    δ = 10^(-8) # avoid comparison of small numbers
    
    # ----- 3d simple cubic ----
    sc3d = gen_kGrid("3dsc-1.2",4)
    ϵk = expandKArr(sc3d, dispersion(sc3d))[:]
    ϵk[abs.(ϵk) .< δ] .= zero(eltype(ϵk))
    
    ϵkq = ϵ_k_plus_q(sc3d, (0.0,0.0,0.0))
    ϵkq[abs.(ϵkq) .< δ] .= zero(eltype(ϵkq))
    @test all(ϵk .≈ ϵkq)
    
    ϵkq = ϵ_k_plus_q(sc3d, (π, π, π))
    ϵkq[abs.(ϵkq) .< δ] .= zero(eltype(ϵkq))
    @test all(ϵk .≈ -ϵkq)
    
    ϵkq = ϵ_k_plus_q(sc3d, (2*π,2*π,2*π))
    ϵkq[abs.(ϵkq) .< δ] .= zero(eltype(ϵkq))
    @test all(ϵk .≈ ϵkq)

    # ----- 2d simple cubic ----
    sc2d = gen_kGrid("2dsc-1.2",4)
    ϵk = expandKArr(sc2d, dispersion(sc2d))[:]
    ϵk[abs.(ϵk) .< δ] .= zero(eltype(ϵk))
    
    ϵkq = ϵ_k_plus_q(sc2d, (0.0,0.0))
    ϵkq[abs.(ϵkq) .< δ] .= zero(eltype(ϵkq))
    @test all(ϵk .≈ ϵkq)

    ϵkq = ϵ_k_plus_q(sc2d, (π, π))
    ϵkq[abs.(ϵkq) .< δ] .= zero(eltype(ϵkq))
    @test all(ϵk .≈ -ϵkq)
    
    ϵkq = ϵ_k_plus_q(sc2d, (2*π,2*π))
    ϵkq[abs.(ϵkq) .< δ] .= zero(eltype(ϵkq))
    @test all(ϵk .≈ ϵkq)

end