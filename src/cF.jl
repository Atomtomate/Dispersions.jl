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

function conv_Indices(::Type{cF}, D::Int, Ns::Int)
    k0 = Tuple(repeat([0],D))
    m1 = floor.(Int, Tuple(repeat([Ns],D)) ./ 2)
    return k0, m1
end

function reduce_KGrid(::Type{cF}, D::Int, Ns::Int, kGrid::AbstractArray)
    (D != 3) && throw(
        ArgumentError(
            "FCC lattice only exists in 3 dimensions!",
        ),
    )
    fsymm(kInd) = fccSymmetries(kInd,Ns)
    ind = collect(Base.product([1:Ns for Di = 1:3]...))
    parents, ops = find_classes(fsymm, vec(ind), UInt32.(repeat([1],48)));
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
    k0, m1 = conv_Indices(cF, D, Ns)
    ind_red = CartesianIndex.(ind_red)
    ind_red_conv  = CartesianIndex.(circshift(ind, m1)[ind_red]); # indices after conv
    ind_red_crossc = CartesianIndex.(circshift(reverse(ind), k0)[ind_red]); # indices after crossc
    # Change from CartesianIndices to LinearIndices for performance reasons
    I = LinearIndices(ind)
    ind_red = I[ind_red]
    ind_red_conv = I[ind_red_conv]
    ind_red_crossc = I[ind_red_crossc]
    expand_perms = map(x -> I[x], expand_perms)
    return ind_red, ind_red_conv, ind_red_crossc, kMult, expand_perms, grid_red
end


function gen_ϵkGrid(::Type{cF}, kGrid::GridPoints, t::T, tp::T, tpp::T) where {T<:Real}
    if tp != 0.0 || tpp != 0.0
        throw(ArgumentError("Dispersion of cF not implemented for non-zero next nearest neightbor hopping!"))
    end
    gen_ϵkGrid(cF, kGrid, t) 
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
        perms = unique(permutations(Tuple(redInd)))
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

function fccSymmetries(kind,Ns)    
    symm = Array{NTuple{3,Int64},1}(undef, 48)
    perms = collect(permutations(kind .-1))
    for i in 1:length(perms)
        symm[8*(i-1)+1] = Tuple(perms[i] .+1)
        symm[8*(i-1)+2] = Tuple((mod.((-perms[i][2]+perms[i][3],-perms[i][1]+perms[i][3],perms[i][3]),Ns)) .+1)
        symm[8*(i-1)+3] = Tuple((mod.((perms[i][1],perms[i][1]-perms[i][3],perms[i][1]-perms[i][2]),Ns)) .+1)
        symm[8*(i-1)+4] = Tuple((mod.((perms[i][2]-perms[i][3],perms[i][2],-perms[i][1]+perms[i][2]),Ns)) .+1)
        symm[8*(i-1)+5] = Tuple((mod.((-perms[i][2]+perms[i][3],-perms[i][2],perms[i][1]-perms[i][2]),Ns)) .+1)
        symm[8*(i-1)+6] = Tuple((mod.((-perms[i][1],-perms[i][1]+perms[i][3],-perms[i][1]+perms[i][2]),Ns)) .+1)
        symm[8*(i-1)+7] = Tuple((mod.((perms[i][2]-perms[i][3],perms[i][1]-perms[i][3],-perms[i][3]),Ns)) .+1)
        symm[8*(i-1)+8] = Tuple((mod.((-perms[i][1],-perms[i][2],-perms[i][3]),Ns)) .+1)
    end
    return unique(symm)[:]
end
