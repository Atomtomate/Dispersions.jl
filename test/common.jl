using FFTW
include("helper_functions.jl")

@testset "reduce_expand" begin
    for gr in grid_list
        for NN = 2:2:8
            kG = gen_kGrid(gr, 4)
            arr0 = randn(rng, ComplexF64, gridshape(kG))
            arr1 = reduceKArr(kG, arr0)
            arr1_exp = expandKArr(kG, arr1)
            expandKArr!(kG, arr1)
            @test all(kG.cache1 .≈ arr1_exp)
            @test all(
                expandKArr(kG, arr1) .≈
                expandKArr(kG, reduceKArr(kG, expandKArr(kG, arr1))),
            ) # consistent?
            if (gr == "3Dsc-1.3"||gr == "fcc-1.4")
                @test sum(kG.kMult) == 4^3
            end
            #TODO: test expand disp and direct disp calc
        end
    end
end

@testset "common functions" begin
    for Ns in [4]
    for gr in grid_list
        kG = gen_kGrid(gr, Ns)
        arr1 = randn(rng, ComplexF64, gridshape(kG))
        arr2 = randn(rng, ComplexF64, gridshape(kG))
        arr1_sym = expandKArr(kG, reduceKArr(kG, arr1))
        arr2_sym = expandKArr(kG, reduceKArr(kG, arr2))
        arr1_fft = fft(arr1)
        arr1_r_fft = fft(reverse(arr1))
        arr2_r_fft = fft(reverse(arr2))
        arr2_fft = fft(arr2)
        conv_theo_res = ifft(arr1_fft .* arr2_fft)
        conv_naive_res,_,_ = naive_conv(arr1, arr2, kG.k0)                     # ∑_k arr1[k] * arr2[q+k]
        conv_theo_res_3 = reverse(ifft(arr1_r_fft .* arr2_fft))
        t1,_,_ = naive_conv(arr1_sym, arr2_sym, kG.k0)
        conv_res = conv(kG, reduceKArr(kG, arr1_sym), reduceKArr(kG, arr2_sym))
        conv_noPlan_res = conv_noPlan(kG, reduceKArr(kG, arr1_sym), reduceKArr(kG, arr2_sym))

        v_full = collect(Dispersions.gen_sampling(grid_type(kG), grid_dimension(kG),  kG.Ns))
        _, dbg_i_pp, dbg_val_pp = naive_conv(v_full, v_full, kG.k0, pp=true, round_entries=false);
        _, dbg_i, dbg_val = naive_conv(v_full, v_full, kG.k0, pp=false, round_entries=false);
        @testset "$gr" begin
            @testset "convolutions" begin
                conv_ph_naive,_,_ = naive_conv(arr1_sym, arr2_sym, kG.k0, pp=false);
                conv_test = Dispersions.conv(kG, reduceKArr(kG,arr1_sym), reduceKArr(kG,arr2_sym), crosscorrelation=true)

                @test all(conv_ph_naive .≈ expandKArr(kG, conv_test) .* Nk(kG))

                conv_pp_naive,_,_ = naive_conv(arr1_sym, arr2_sym, kG.k0, pp=true);
                conv_pp_test = Dispersions.conv(kG, reduceKArr(kG,arr1_sym), reduceKArr(kG, arr2_sym), crosscorrelation=false)
                @test all(conv_pp_naive ./ kG.Nk .≈ expandKArr(kG, conv_pp_test))


                @test all(reduceKArr(kG, expandKArr(kG, conv_res)) .≈ conv_res)     # is symmetry preserved?

                q_list = map(el->check_q(el, Nk(kG)), dbg_val);
                # preprocess rounds all entries in order for unqiue to work
                q_list_pre = map(x -> map(xi -> round.(xi,digits=8), x), q_list)
                q_list_pre = map(x -> map(xi -> map(xii -> xii ≈ 0 ? abs(xii) : xii, xi), x), q_list_pre)
                # Test if all entries have the same common q difference
                @test all(map(qi_list -> length(unique(qi_list)) == 1, q_list_pre))
                # Extract only q vector (since we checked that only one unique exists
                q_list_single = map(x->x[1], q_list_pre)
                
                # bcc and fcc are only basis transformed simple cubic lattices
                if !any(map(s -> contains(gr,s), ["bcc", "fcc", "cI", "cF"]))
                    # Test if all q-vectors are at the correct index
                    @test all(map(x-> all(x[1] .≈ x[2]), zip(q_list_single, v_full)))
                    # Test if all possible q-vectors are contained in convolution
                    @test all(map(el -> all(el[1] .≈ el[2]), zip(sort(unique(v_full)[:]), sort(map(qi->qi[1], q_list_pre)[:]))))
                    # Test if all combinations of k_i and k_{i+j} are contained in each result index
                    @test all(map(entry -> length(unique(map(x->x[1], entry))), dbg_val) .== Nk(kG))
                    @test all(map(entry -> length(unique(map(x->x[2], entry))), dbg_val) .== Nk(kG))

                    conv_ph_fft_test = circshift(reverse(ifft(fft(arr1_sym) .* fft(reverse(arr2_sym)))), kG.k0 .- 1)
                    @test all(conv_ph_fft_test .≈ conv_ph_naive)
                end

                q_list_pp = map(el->check_q(el, Nk(kG), pp=true), dbg_val_pp);
                # preprocess rounds all entries in order for unqiue to work
                q_list_pp_pre = map(x -> map(xi -> round.(xi,digits=8), x), q_list_pp)
                q_list_pp_pre = map(x -> map(xi -> map(xii -> xii ≈ 0 ? abs(xii) : xii, xi), x), q_list_pp_pre)
                # Test if all entries have the same common q difference
                @test all(map(qi_list -> length(unique(qi_list)) == 1, q_list_pp_pre))
                # extract only q vector (since we checked that only one unique exists
                q_list_pp_single = map(x->x[1], q_list_pp_pre)

                # bcc and fcc are only basis transformed simple cubic lattices
                if !any(map(s -> contains(gr,s), ["bcc", "fcc", "cI", "cF"]))
                    # test if all q-vectors are at the correct index
                    @test all(map(x-> all(x[1] .≈ x[2]), zip(q_list_pp_single, v_full)))
                    # Test if all possible q-vectors are contained in convolution
                    @test all(map(el -> all(el[1] .≈ el[2]), zip(sort(unique(v_full)[:]), sort(map(qi->qi[1], q_list_pp_pre)[:]))))
                    # Test if all combinations of k_i and k_{i+j} are contained in each result index
                    @test all(map(entry -> length(unique(map(x->x[1], entry))), dbg_val_pp) .== Nk(kG))
                    @test all(map(entry -> length(unique(map(x->x[2], entry))), dbg_val_pp) .== Nk(kG))

                    conv_pp_fft_test = circshift(ifft(fft(arr1_sym) .* fft(arr2_sym)), -1 .* kG.k0 .+ 1)
                    @test all(conv_pp_fft_test .≈ conv_pp_naive)
                end
            end
            @testset "reverse" begin
                # bcc and fcc are only basis transformed simple cubic lattices
                if !any(map(s -> contains(gr,s), ["bcc", "fcc", "cI", "cF"]))
                    minus_vec_naive = map(k -> transform_k_to_minus_k(k, kG), v_full)
                    minus_v_full = reverseKArr(kG, v_full)
                    #minus_v_full_red = reduceKArr(kG, minus_v_full)
                    @test all(map(x-> all(x[1] .≈ x[2]), zip(minus_vec_naive, minus_v_full)))
                    
                    #TODO: this needs more tests and impl. for other lattices!
                    zero_vec = Tuple(repeat([0],grid_dimension(kG)))
                    @test all(transform_to_first_BZ(kG, zero_vec) .≈ zero_vec)
                end
            end
        end
        end
    end
end

@testset "symmetry paths" begin
end

