@testset "helpers" begin
    for Ns in [4,6]
        for gr in grid_list
            kG = gen_kGrid(gr, Ns)

            @testset "build k minus kp" begin
                #TODO: not implemented
                if !any(map(s -> contains(gr,s), ["bcc", "fcc", "cI", "cF"]))
                v_full = gen_sampling(grid_type(kG), grid_dimension(kG), kG.Ns)
                q_lookup = build_q_lookup(kG)
                fails = []
                for k in v_full
                    for kp in v_full
                        q_test = round.(transform_to_first_BZ(kG, k .- kp), digits=6)
                        !(q_test in keys(q_lookup)) && push!(fails, 1)
                    end
                end
                !isempty(fails) && println(gr, " // ", Ns)
                @test isempty(fails)
                end
            end
        end
    end
end

