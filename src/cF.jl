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
    (D != 3) && throw(
        ArgumentError(
            "FCC lattice only exists in 3 dimensions!",
        ),
    )
    ind = collect(Base.product([1:Ns for Di = 1:D]...))
    parents, ops = find_classes(fccSymmetries, kGrid[:], UInt32.(repeat([1],48)));
    kmap, grid_red = minimal_set(parents, kGrid[:]);  
    kMult = zeros(Int64,length(grid_red))
    expand_perms = Vector{Vector{NTuple{D, Int}}}(undef,length(grid_red))
    for i in 1:length(grid_red)
        kMult[i] = length(fccSymmetries(grid_red[i]))
        expand_perms[i] = fccSymmetries(grid_red[i])
    end

    #ind = collect(Base.product([1:Ns for Di = 1:3]...))
    #kMult = ones(length(ind))
    #expand_perms = map(x -> [CartesianIndex{3}(x)], ind[:])
    red_map = CartesianIndex.(ind[:])
    red_conv_map = reverse(red_map)
    return red_map, red_conv_map, kMult, expand_perms, grid_red
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

function fccSymmetries(k)
	symm = Array{Vector{Float64},1}(undef, 48)
	perms = collect(permutations(k))
	for i in 1:length(perms)
		symm[8*(i-1)+1] = perms[i] .* [1,1,1]
		symm[8*(i-1)+2] = perms[i] .* [1,1,-1]
		symm[8*(i-1)+3] = perms[i] .* [1,-1,1]
		symm[8*(i-1)+4] = perms[i] .* [1,-1,-1]
		symm[8*(i-1)+5] = perms[i] .* [-1,1,1]
		symm[8*(i-1)+6] = perms[i] .* [-1,1,-1]
		symm[8*(i-1)+7] = perms[i] .* [-1,-1,1]
		symm[8*(i-1)+8] = perms[i] .* [-1,-1,-1]
	end
	return unique(symm)[:]
end