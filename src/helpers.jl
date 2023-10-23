
"""
    build_q_lookup(kG::KGrid)

Builds lookup table of indices, for all q-vectors in `kG`.
"""
function build_q_lookup(kG::KGrid)

    v_full = gen_sampling(grid_type(kG), grid_dimension(kG), kG.Ns)
    v_full_rounded = map(ki -> round.(ki, digits=6), v_full)
    v_full_binned  = map(li -> map(ind -> v_full_rounded[ind], li), kG.expand_perms)

    q_lookup = Dict{eltype(v_full_rounded), Int}()
    sizehint!(q_lookup, length(v_full_rounded))

    for (ki,red_q_list) in enumerate(v_full_binned)
        for qi in red_q_list
            q_lookup[qi] = ki
        end
    end

    return q_lookup
end
