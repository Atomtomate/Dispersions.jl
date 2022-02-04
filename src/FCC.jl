#TODO update tests and docu for improved expandKArr


"""
    gen_cF_kGrid(Nk::Int64, t::Float64, sampling::AbstractArray)

Generates a face centered cubic lattice of type k-grid
See also [`FullKGrid_cF_3D`](@ref)
# Examples
```
julia> gen_cF_kGrid(2, 2)
"""
function gen_cF_kGrid(Nk::Int64, t::Float64, sampling::AbstractArray)
	return FullKGrid_cF(Nk,t,sampling)
end

# ================================================================================ #
#                               Face Centered Cubic                                #
# ================================================================================ #

# -------------------------------------------------------------------------------- #
#                                     Types                                        #
# -------------------------------------------------------------------------------- #

abstract type cF <: KGridType end

"""
    FullKGrid_cF  <: FullKGrid{cF}

Fields
-------------
- **`kGrid`** : `Array{Tuple{Float64, ...}}` of kGrids. Each element is a 3-tuple
"""
struct FullKGrid_cF <: FullKGrid{cF, 3}
    Nk::Int
    Ns::Int
    kGrid::Array{NTuple{3,Float64},1}
    ϵkGrid::GridDisp
    t::Float64
    fftw_plan::FFTW.cFFTWPlan
    function FullKGrid_cF(Nk::Int, t::Float64, sampling::AbstractArray; fftw_plan=nothing)
        kGrid  = collect(Base.product([sampling for Di in 1:3]...))[:]
        fftw_plan = fftw_plan === nothing ? plan_fft!(randn(Complex{Float64}, repeat([Nk], 3)...), flags=FFTW.ESTIMATE, timelimit=Inf) : fftw_plan
        new(Nk^3, Nk, kGrid, gen_ϵkGrid(cF,kGrid,t),t,fftw_plan)
    end
end

"""
    ReducedKGrid_cF  <: ReducedKGrid{cF}

Reduced k grid only containing points necessary for computation of quantities involving 
this grid type and multiplicity of each stored point.

Fields
-------------
- **`Ns`**     		 : `Int` k points in each dimension.
- **`Nk`**     		 : `Int` Total number of k points.
- **`kInd`**   		 : `Array{Tuple{Int,...}}` indices in full grid. used for reconstruction of full grid.
- **`kMult`**  		 : `Array{Float64,1}` multiplicity of point, used for calculations involving reduced k grids.
- **`kGrid`**  		 : `Array{Tuple{Float64,...}}` k points of reduced grid.
- **`ϵkGrid`** 		 : `GridDisp` Dispersion at each k point.
- **`expand_perms`** : List of length `kMult` lists. List `i` contains all indices, to which the i-th point in the reduced grid is mapped, when exapnding to the full grid.
- **`expand_cache`** : `Array{Complex{Float64}}` Cache for performance reasons.
- **`fftw_plan`** 	 : `FFTW.cFFTWPlan` Optional plan for the FFTW.
"""
struct ReducedKGrid_cF  <: ReducedKGrid{cF,3}
    Nk::Int
    Ns::Int
    t::Float64
    kInd::Array{NTuple{3,Int},1}
    kMult::Array{Float64,1}
    kGrid::Array{NTuple{3,Float64},1}
    ϵkGrid::GridDisp
    expand_perms::Vector{Vector{CartesianIndex{3}}}
    expand_cache::Array{Complex{Float64}}
    fftw_plan::FFTW.cFFTWPlan
end

# -------------------------------------------------------------------------------- #
#                                   Interface                                      #
# -------------------------------------------------------------------------------- #

gridshape(kG::ReducedKGrid_cF) = ntuple(_ -> kG.Ns, 3)
gridshape(kG::FullKGrid_cF) = ntuple(_ -> kG.Ns, 3)

# ---------------------------- BZ to f. irr. BZ -------------------------------

"""
    reduceKGrid(kG::FullKGrid{T}) where T <: cF

Returns the grid on the fully irredrucible BZ.
Filters an arbitrary grid, defined on a full kGrid, so that only 
the lower triangle remains, i.e.
for any (x_1, x_2, ...) the condition x_1 >= x_2 >= x_3 ... 
is fulfilled.
"""
function reduceKGrid(kG::FullKGrid{cF,3})
    kGrid = reshape(kG.kGrid, gridshape(kG))
    ϵkGrid = reshape(kG.ϵkGrid, gridshape(kG))
    ind = collect(Base.product([1:kG.Ns for Di in 1:3]...))
	D = 3
    ll = floor(Int,size(kGrid,1)/2 + 1) - 1
    la = ceil(Int,kG.Ns/2)
    index = if D == 2
        [[x,y]  for x=la:ll+la for y=la:x]
    elseif D == 3
        [[x,y,z] for x=la:ll+la for y=la:x for z = la:y]
    else
        error("D ∉ [2,3] not implemented yet!")
    end

    ind_red = GridInd{3}(undef, length(index))
    grid_red = GridPoints{3}(undef, length(index))
    ϵk_red = GridDisp(undef, length(index))

    # Compute reduced arrays
    for (i,ti) in enumerate(index)
        ind_red[i] = ind[ti...]
        grid_red[i] = kGrid[ti...]
        ϵk_red[i] = ϵkGrid[ti...]
    end

    kMult, expand_perms = build_expand_mapping_SC(D, kG.Ns, ind_red)
    expand_cache = Array{Complex{Float64}}(undef, gridshape(kG))

    return ReducedKGrid_cF(kG.Nk, kG.Ns, kG.t, ind_red, kMult, grid_red, ϵk_red,
                              expand_perms, expand_cache, kG.fftw_plan)
end

"""
    expandKArr(kG::ReducedKGrid{T1}, arr::Array{T2,1})

Expands array of values on reduced k grid back to full BZ.
"""
#TODO: refactor, these functions are defined almost equivalently in SC,hex
function expandKArr(kG::ReducedKGrid_cF, arr::AbstractArray{T, 1})::AbstractArray{T, 3} where {T}
    length(arr) != length(kG.kInd) && throw(ArgumentError("length of k grid ($(length(kG.kInd))) and argument ($(length(arr))) not matching"))
    res = similar(arr, gridshape(kG)...)
    expandKArr!(kG, res, arr)
    return res
end

function expandKArr!(kG::ReducedKGrid_cF, arr::Array{Complex{Float64}, 1})
    for (ri,perms) in enumerate(kG.expand_perms)
        @simd for p in perms
            @inbounds kG.expand_cache[p] = arr[ri]
        end
    end
end

function expandKArr!(kG::ReducedKGrid_cF, res::AbstractArray{T,3}, arr::Array{T, 1}) where {T}
    for (ri,perms) in enumerate(kG.expand_perms)
        @simd for p in perms
            @inbounds res[p] = arr[ri]
        end
    end
end

function reduceKArr(kG::ReducedKGrid_cF, arr::AbstractArray{T,3}) where {T}
    res = similar(arr, length(kG.kInd))
    reduceKArr!(kG, res, arr)
    return res
end

function reduceKArr!(kG::ReducedKGrid_cF, res::AbstractArray{T,1}, arr::AbstractArray{T,3}) where {T}
    for (i,ki) in enumerate(kG.kInd)
        @inbounds res[i] = arr[ki...]
    end
end

#TODO: implement
gen_ϵkGrid(::Type{cF}, kGrid::GridPoints, t::T) where T <: Real = collect(map(kᵢ -> -2*t*sum(cos.(kᵢ)), kGrid))
ifft_post(kG::ReducedKGrid_cF, x::Array{T,N}) where {N, T <: Number} = ShiftedArrays.circshift(x, floor.(Int, gridshape(kG) ./ 2) .+ 1)

function build_expand_mapping_cF(D::Int, Ns::Int, ind_red::Array)
    expand_perms = Vector{Vector{CartesianIndex{3}}}(undef, length(ind_red))
    kMult = Array{Int,1}(undef, length(ind_red))
	D = 3
    mirror_list = Array{NTuple{D,Int},1}()
    for i in 1:D
        push!(mirror_list, map( x-> tuple((x)...) ,unique(permutations([(j <= i) ? Ns : 0 for j in 1:D ])))...)
    end
    #  - Expand mapping
    for (ri, redInd) in enumerate(ind_red)
        perms = unique(permutations(redInd))
        expand_perms[ri] = Vector{CartesianIndex{D}}()
        for (ip,p) in enumerate(perms)
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
