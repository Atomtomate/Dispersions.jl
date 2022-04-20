# ================================================================================ #
#                               Face Centered Cubic                                #
# ================================================================================ #

abstract type cF <: KGridType end

# -------------------------------------------------------------------------------- #
#                            Grid Generator Functions                              #
# -------------------------------------------------------------------------------- #

gen_sampling(::Type{cF}, D::Int, Ns::Int) =
    Base.product([[(2 * π / Ns) * (j - 1) for j = 1:Ns] for Di = 1:3]...)
basis_transform(::Type{cF}, v::Tuple) =
    Tuple([-1.0 1.0 1.0; 1.0 -1.0 1.0; 1.0 1.0 -1.0] * collect(v))

function reduce_KGrid(::Type{cF}, D::Int, Ns::Int, kGrid::AbstractArray)
    ind = collect(Base.product([1:Ns for Di = 1:3]...))
    kMult = ones(length(ind))
    expand_perms = map(x -> [CartesianIndex{3}(x)], ind[:])
    return CartesianIndex.(ind[:]), kMult, expand_perms, kGrid[:]
end

gen_ϵkGrid(::Type{cF}, kGrid::GridPoints, t::T) where {T<:Real} = collect(
    map(
        kᵢ ->
            -4t*
            (cos(kᵢ[1]/2) * cos(kᵢ[2]/2) + cos(kᵢ[1]/2) * cos(kᵢ[3]/2) + cos(kᵢ[2]/2) * cos(kᵢ[3]/2)),
        kGrid,
    ),
)

# -------------------------------------------------------------------------------- #
#                             Custom Helper Functions                              #
# -------------------------------------------------------------------------------- #

conv_sample_post(kG::KGrid{cF,3}, x) = x

function build_expand_mapping_cF(D::Int, Ns::Int, ind_red::Array)
    expand_perms = Vector{Vector{CartesianIndex{3}}}(undef, length(ind_red))
    kMult = Array{Int,1}(undef, length(ind_red))
    D = 3
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
