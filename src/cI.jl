# ================================================================================ #
#                               Body Centered Cubic                                #
# ================================================================================ #

abstract type cI <: KGridType end

# -------------------------------------------------------------------------------- #
#                            Grid Generator Functions                              #
# -------------------------------------------------------------------------------- #

gen_sampling(::Type{cI}, D::Int, Ns::Int) =
    Base.product([[(2 * π / Ns) * (j - 1) for j = 1:Ns] for Di = 1:3]...)
basis_transform(::Type{cI}, v::Tuple) =
    Tuple([1.0 1.0 0.0; 0.0 1.0 1.0; 1.0 0.0 1.0] * collect(v))

function reduce_KGrid(::Type{cI}, D::Int, Ns::Int, kGrid::AbstractArray)
    (D != 3) && throw(
        ArgumentError(
            "BCC lattice only exists in 3 dimensions!",
        ),
    )
    fsymm(kInd) = bccSymmetries(kInd)
    ind = collect(Base.product([1:Ns for Di = 1:3]...))
    parents, ops = find_classes(fsymm, vec(ind), UInt32.(repeat([1],12)));
    kmap, ind_red = minimal_set(parents, vec(ind))
    grid_red = Array{NTuple{3,Float64},1}(undef,length(ind_red))
    for (i,indi) in enumerate(ind_red)
        grid_red[i] = kGrid[CartesianIndex(indi)]
    end
    kMult = zeros(Int64,length(ind_red))
    expand_perms = Vector{Vector{CartesianIndex{D}}}(undef,length(ind_red))
    for i in 1:length(ind_red)
        kMult[i] = length(fsymm(ind_red[i]))
        expand_perms[i] = CartesianIndex.(fsymm(ind_red[i]))
    end
    red_conv_map = CartesianIndex.(reverse(ind)[CartesianIndex.(ind_red)])
    return CartesianIndex.(ind_red), red_conv_map, kMult, expand_perms, grid_red
end

function gen_ϵkGrid(::Type{cI}, kGrid::GridPoints, t::T, tp::T, tpp::T) where {T<:Real}
    if tp != 0.0 || tpp != 0.0
        throw(ArgumentError("Dispersion of cF not implemented for non-zero next nearest neightbor hopping!"))
    end
    gen_ϵkGrid(cF, kGrid, t) 
end

gen_ϵkGrid(::Type{cI}, kGrid::GridPoints, t::T) where {T<:Real} = collect(
    map(
        kᵢ ->
            -8t*
            (cos(kᵢ[1]/2) * cos(kᵢ[2]/2) * cos(kᵢ[3]/2)),
        kGrid,
    ),
)

# -------------------------------------------------------------------------------- #
#                             Custom Helper Functions                              #
# -------------------------------------------------------------------------------- #

conv_sample_post(kG::KGrid{cI,3}, x) = x

function build_expand_mapping_cI(D::Int, Ns::Int, ind_red::Array)
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

function bccSymmetries(kind)
    return unique(Tuple.(collect(permutations(kind))))[:]
end
