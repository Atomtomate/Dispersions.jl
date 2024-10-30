# ================================================================================ #
#                                   Simple Cubic                                   #
# ================================================================================ #

abstract type cP <: KGridType end

# -------------------------------------------------------------------------------- #
#                            Grid Generator Functions                              #
# -------------------------------------------------------------------------------- #

gen_sampling(::Type{cP}, D::Int, Ns::Int) =
    Base.product([[(2 * π / Ns) * j - π for j = 1:Ns] for Di = 1:D]...)

basis_transform(::Type{cP}, v::Tuple) = v

transform_to_first_BZ(kG::KGrid{cP,D}, k) where D =
    Tuple(map(ki -> ki != 0 ? mod(ki + (π - 1/kG.Nk), 2π) + (-π) + 1/kG.Nk : ki, k))
    
function conv_Indices(::Type{cP}, D::Int, Ns::Int)
    k0 = floor.(Int, Tuple(repeat([Ns],D)) ./ 2) .- 1
    m1 = -1 .* k0 .+ 0
    return k0, m1
end

function reduce_KGrid(::Type{cP}, D::Int, Ns::Int, kGrid::AbstractArray)
    (Ns % 2 != 0 && Ns != 1) && 
        println("WARNING! k-sampling with uneven grid-points may break functionality!")
        
    ind = collect(Base.product([1:Ns for Di = 1:D]...))
    ll = floor(Int, size(kGrid, 1) / 2 + 1) - 1
    la = ceil(Int, Ns / 2)
    index = if D == 2
        [CartesianIndex(x, y) for x = la:ll+la for y = la:x]
    elseif D == 3
        [CartesianIndex(x, y, z) for x = la:ll+la for y = la:x for z = la:y]
    elseif D == 4
        [CartesianIndex(x1, x2, x3, x4) for x1 = la:ll+la for x2 = la:x1 for x3 = la:x2 for x4 = la:x3]
    else
        error("cP for D ∉ [2,3,4] not implemented yet!")
    end

    ind_red = Vector{CartesianIndex{D}}(undef, length(index))
    grid_red = GridPoints{D}(undef, length(index))

    # Compute reduced arrays
    for (i, ti) in enumerate(index)
        ind_red[i] = CartesianIndex(ind[ti])
        grid_red[i] = kGrid[ti]
    end

    k0, m1 = conv_Indices(cP, D, Ns)

    ind_red_conv  = CartesianIndex.(circshift(ind, m1)[index]); # indices after conv
    ind_red_crossc = CartesianIndex.(circshift(reverse(ind), k0)[index]); # indices after crossc
    kMult, expand_perms = build_expand_mapping_cP(D, Ns, ind_red)

    # Change from CartesianIndices to LinearIndices for performance reasons
    I = LinearIndices(ind)
    index = I[index]
    ind_red_conv = I[ind_red_conv]
    ind_red_crossc = I[ind_red_crossc]
    return index, ind_red_conv, ind_red_crossc, kMult, expand_perms, grid_red
end

function gen_ϵkGrid(::Type{cP}, kGrid::GridPoints, t::T, tp::T, tpp::T) where {T<:Real}
    if tp != 0.0 || tpp != 0.0
        throw(ArgumentError("Dispersion of cP not implemented for non-zero next nearest neightbor hopping!"))
    end
    gen_ϵkGrid(cP, kGrid, t) 
end

gen_ϵkGrid(::Type{cP}, kGrid::GridPoints, t::T) where {T<:Real} =
    collect(map(kᵢ -> -2 * t * sum(cos.(kᵢ)), kGrid))

# -------------------------------------------------------------------------------- #
#                             Custom Helper Functions                              #
# -------------------------------------------------------------------------------- #
function symmetry_points(kG::KGrid)
    [()]
end

"""
    conv_sample_post(kG::KGrid{cP,D}, x)

This is needed in case the sampling is not starting at 0, i.e. [0,V) × [0,V) ... in order to shift the 0 frequency to the appropriate sampling point.   
"""
conv_sample_post(kG::KGrid{cP,D}, x; crosscorrelation=true) where {D} =
    ShiftedArrays.circshift(x, floor.(Int, gridshape(kG) ./ 2) .- 1)
#TODO: this somehow works when not doing the reverse on the second input. We should find out why, this makes the convolution a lot faster
conv_post_old(kG::KGrid{cP,D}, x::Array{T,D}) where {D,T<:Number} =
    reduceKArr(kG, ShiftedArrays.circshift(x, floor.(Int, gridshape(kG) ./ 2) .+ 1))

function build_expand_mapping_cP(D::Int, Ns::Int, ind_red::Array)
    expand_perms = Vector{Vector{CartesianIndex{D}}}(undef, length(ind_red))
    kMult = Array{Int,1}(undef, length(ind_red))

    mirror_list = Array{NTuple{D,Int},1}()
   for i = 1:D
        push!(
            mirror_list,
            map(
                x -> tuple((x)...),
                unique(permutations([(j <= i) ? Ns : 0 for j = 1:D])),
            )...,
        )
    end
    #  - Expand mapping
    for (ri, redInd) in enumerate(ind_red)
        perms = unique(permutations(redInd))
        expand_perms[ri] = Vector{CartesianIndex{D}}()
        for (ip, p) in enumerate(perms)
            push!(expand_perms[ri], CartesianIndex(p...))
            for mi in mirror_list
                if all(abs.(mi .- p) .> 0)
                    push!(expand_perms[ri], CartesianIndex(abs.(mi .- p)...))
                end
            end
        end
        expand_perms[ri] = unique(expand_perms[ri])
        kMult[ri] = length(expand_perms[ri])
    end
    return kMult, expand_perms
end
