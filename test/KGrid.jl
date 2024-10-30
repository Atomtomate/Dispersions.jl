using Base.Iterators
# TODO: test kvices/ifft
# TODO: test ϵ_k_plus_q for other grid types than simple cubic

@testset "gen_kGrid" begin
    for NN in [4, 6, 8]
        gl = map(x -> gen_kGrid(x, NN), grid_list)
        for i = 1:length(gl)
            kG = gl[i]
            Di = grid_dimension(kG)
            @test kG.Nk == NN^Di
            @test kG.Ns == NN
            @test all(size(kG.cache1) .== repeat([NN], Di))
            @test all(size(kG.cache2) .== repeat([NN], Di)) 
        end
    end
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

@testset "grid properties" begin
    # cF
    fcc = gen_kGrid("fcc-1.2",2)
    @test grid_dimension(fcc) == 3
    @test grid_type(fcc) === cF
    
    # cI
    bcc = gen_kGrid("bcc-1.2",2)
    @test grid_dimension(bcc) == 3
    @test grid_type(bcc) === cI

    # cP
    sc2d = gen_kGrid("2Dsc-1.3",2)
    @test grid_dimension(sc2d) == 2
    @test grid_type(sc2d) === cP
    
    sc3d = gen_kGrid("3Dsc-1.3",2)
    @test grid_type(sc3d) === cP
    @test grid_dimension(sc3d) == 3

    # cPnn
    sc2dnn = gen_kGrid("2Dsc-1.3--1.4-1.5",2)
    @test grid_dimension(sc2dnn) == 2
    @test grid_type(sc2dnn) === cPnn
    
end

@testset "grid subsampling" begin
    #TODO: this is only a first draft implementation
    kG_full = gen_kGrid("2Dsc-1.1", 20)
    kG_sub,ii = build_kGrid_subsample(kG_full, 11)
    tt = map(x-> round.(x,digits=4), dispersion(kG_full))
    @test all(map(x-> round.(x,digits=4) in tt, dispersion(kG_sub)))
    @test kG_sub.Ns == 10
end
