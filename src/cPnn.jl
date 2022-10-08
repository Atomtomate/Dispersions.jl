# ================================================================================ #
#                                   Simple Cubic                                   #
# ================================================================================ #

abstract type cPnn <: KGridType end

# -------------------------------------------------------------------------------- #
#                            Grid Generator Functions                              #
# -------------------------------------------------------------------------------- #

gen_sampling(::Type{cPnn}, D::Int, Ns::Int) =
    Base.product([[(2 * π / Ns) * j - π for j = 1:Ns] for Di = 1:D]...)
basis_transform(::Type{cPnn}, v::Tuple) = v

function reduce_KGrid(::Type{cPnn}, D::Int, Ns::Int, kGrid::AbstractArray)
    ind = collect(Base.product([1:Ns for Di in 1:2]...))
    kMult = ones(length(ind))
    expand_perms = map(x -> [CartesianIndex{2}(x)],ind[:])
    red_map = CartesianIndex.(ind[:])
    red_conv_map = reverse(red_map)
    return red_map, red_conv_map, kMult, expand_perms, kGrid[:]
end

gen_ϵkGrid(::Type{cPnn}, kGrid::GridPoints, t::T) where {T<:Real} =
    collect(map(kᵢ -> -2 * t * sum(cos.(kᵢ)), kGrid))

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
