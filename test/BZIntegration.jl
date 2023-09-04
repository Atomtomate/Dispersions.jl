@testset "KSum" begin
    for grid in grid_list
        for NN = 2:2:8
            gr = gen_kGrid(grid, NN)
            @test kintegrate(gr, ones(size(gr.ϵkGrid)))[1] ≈ 1.0
            @test_throws ArgumentError kintegrate(gr, ones(7, 18), 2)
            @test all(kintegrate(gr, ones(length(gr.ϵkGrid)..., 4), 1) .≈ ones(4))
        end
    end

    # --- check against analytic properties of the 3d simple cubic lattice ---
    # Integral of the dispersion raised to odd powers vanishes. This can be shown by induction. For the indution step one may use integration by reduction.
    kGr = gen_kGrid("3dsc-$(1/(2*sqrt(6)))", 20)
    @test abs(kintegrate(kGr, kGr.ϵkGrid))    < 10^(-10) 
    @test abs(kintegrate(kGr, kGr.ϵkGrid.^3)) < 10^(-10)
    @test abs(kintegrate(kGr, kGr.ϵkGrid.^5)) < 10^(-10)
    # Integral of the dispersion raised to even powers can be evaluated by hand
    @test kintegrate(kGr, kGr.ϵkGrid.^2) ≈ 1/4
    @test kintegrate(kGr, kGr.ϵkGrid.^4) ≈ 5/32
    @test kintegrate(kGr, kGr.ϵkGrid.^6) ≈ 155/1152
end
