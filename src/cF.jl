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
    fsymm(kInd) = fccSymmetries(kInd,Ns)
    ind = collect(Base.product([1:Ns for Di = 1:3]...))
    parents, ops = find_classes(fsymm, vec(ind), UInt32.(repeat([1],12)));
    kmap, ind_red = minimal_set(parents, vec(ind));
    println(ind_red)
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
    #ind = collect(Base.product([1:Ns for Di = 1:3]...))
    #kMult = ones(length(ind))
    #expand_perms = map(x -> [CartesianIndex{3}(x)], ind[:])
    #red_map = CartesianIndex.(ind[:])
    red_conv_map = CartesianIndex.(reverse(ind_red))
    return CartesianIndex.(ind_red), red_conv_map, kMult, expand_perms, grid_red
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

function fccSymmetries(kind,Ns)    
    symm = Array{NTuple{3,Int64},1}(undef, 12)
    perms = collect(permutations(kind.-1))
    for i in 1:length(perms)
        symm[2*(i-1)+1] = Tuple(perms[i] .+ 1)
        symm[2*(i-1)+2] = Tuple(((Ns-1) .- perms[i]) .+1)
    end
    return unique(symm)[:]
end