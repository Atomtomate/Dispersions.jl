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

function reduce_KGrid(::Type{cP}, D::Int, Ns::Int, kGrid::AbstractArray)
    (Ns % 2 != 0 && Ns != 1) && throw(
        ArgumentError(
            "Cannot reduce simple cubic lattice with uneven number of samplng points!",
        ),
    )
    ind = collect(Base.product([1:Ns for Di = 1:D]...))
    ll = floor(Int, size(kGrid, 1) / 2 + 1) - 1
    la = ceil(Int, Ns / 2)
    index = if D == 2
        [CartesianIndex(x, y) for x = la:ll+la for y = la:x]
    elseif D == 3
        [CartesianIndex(x, y, z) for x = la:ll+la for y = la:x for z = la:y]
    else
        error("cP for D ∉ [2,3] not implemented yet!")
    end

    ind_red = GridInd{D}(undef, length(index))
    grid_red = GridPoints{D}(undef, length(index))

    # Compute reduced arrays
    for (i, ti) in enumerate(index)
        ind_red[i] = CartesianIndex(ind[ti])
        grid_red[i] = kGrid[ti]
    end

    kMult, expand_perms = build_expand_mapping_cP(D, Ns, ind_red)
    return index, kMult, expand_perms, grid_red
end

gen_ϵkGrid(::Type{cP}, kGrid::GridPoints, t::T) where {T<:Real} =
    collect(map(kᵢ -> -2 * t * sum(cos.(kᵢ)), kGrid))

"""
    conv_post!(kG::KGrid{cP,D}, res::Array{T,1}, x::Array{T,D}) where {D,T} 

Inplace version of [`conv_post`](@ref). Warning: `res` cannot alias kG.cache2!
"""
function conv_post!(kG::KGrid{cP,D}, res::Array{T,1}, x::Array{T,D}) where {D,T} 
    reverse!(x)
    ShiftedArrays.circshift!(kG.cache2, x, floor.(Int, gridshape(kG) ./ 2) .- 1)
    reduceKArr!(kG, res, kG.cache2)
    norm = Nk(kG)
    res[:] = res ./ norm
end
#TODO: optimimize this, i.e. write reduce_from_conv function
#      Description: reverse is part of rewriting the usual to our convolution definition. The circshift rotates the q= 0 point back into the middle of the array (since we sample from -pi to pi for cP.

# -------------------------------------------------------------------------------- #
#                             Custom Helper Functions                              #
# -------------------------------------------------------------------------------- #

conv_sample_post(kG::KGrid{cP,D}, x) where {D} =
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
