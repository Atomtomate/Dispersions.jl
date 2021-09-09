import Base.collect
#TODO update tests and docu for improved expandKArr


"""
    gen_cP_kGrid(Nk::Int64, D::Int64, t::Float64)

Generates a simple cubic lattice in `D` Dimensions
See also [`FullKGrid_cP_3D`](@ref) and [`FullKGrid_cP_3D`](@ref)
# Examples
```
julia> gen_cP_kGrid(2, 2)
"""
function gen_cP_kGrid(Nk::Int64, D::Int64, t::Float64)
    if D == 2
        return FullKGrid_cP(2, Nk, t)
    elseif D == 3
        return FullKGrid_cP(3, Nk, t)
    else
        throw("Simple Cubic only implemented for 2D and 3D")
    end
end


# ================================================================================ #
#                                Simple Cubic 2D                                    #
# ================================================================================ #

# -------------------------------------------------------------------------------- #
#                                     Types                                        #
# -------------------------------------------------------------------------------- #

abstract type cP <: KGridType end

"""
    FullKGrid_cP{D}  <: FullKGrid{cP{D}}

Fields
-------------
- **`kGrid`** : `Array{Tuple{Float64, ...}}` of kGrids. Each element is a D-tuple
"""
struct FullKGrid_cP{D} <: FullKGrid{cP, D}
    Nk::Int
    Ns::Int
    kGrid::Array{NTuple{D,Float64},1}
    ϵkGrid::GridDisp
    t::Float64
    fftw_plan::FFTW.cFFTWPlan
    function FullKGrid_cP(D::Int, Nk::Int, t::Float64; fftw_plan=nothing)
        kx = [(2*π/Nk) * j - π for j in 1:Nk]
        kGrid  = collect(Base.product([kx for Di in 1:D]...))[:]
        fftw_plan = fftw_plan === nothing ? plan_fft!(randn(Complex{Float64}, repeat([Nk], D)...), flags=FFTW.ESTIMATE, timelimit=Inf) : fftw_plan
        new{D}(Nk^D, Nk, kGrid, gen_ϵkGrid(cP,kGrid,t),t,fftw_plan)
    end
end

"""
    ReducedKGrid_cP{D}  <: ReducedKGrid{cP{D}}

Reduced k grid only containing points necessary for computation of quantities involving 
this grid type and multiplicity of each stored point.

Fields
-------------
- **`Nk`**    : `Int` k points in full grid.
- **`kInd`**  : `Array{Tuple{Int,...}}` indices in full grid. used for reconstruction of full grid.
- **`kMult`** : `Array{Float64,1}` multiplicity of point, used for calculations involving reduced k grids.
- **`kGrid`** : `Array{Tuple{Float64,...}}` k points of reduced grid.
"""
struct ReducedKGrid_cP{D}  <: ReducedKGrid{cP,D}
    Nk::Int
    Ns::Int
    t::Float64
    kInd::Array{NTuple{D,Int},1}
    kMult::Array{Float64,1}
    kGrid::Array{NTuple{D,Float64},1}
    ϵkGrid::GridDisp
    expand_perms::Vector{Vector{CartesianIndex{D}}}
    expand_cache::Array{Complex{Float64}}
    fftw_plan::FFTW.cFFTWPlan
end

# -------------------------------------------------------------------------------- #
#                                   Interface                                      #
# -------------------------------------------------------------------------------- #

gridshape(kG::ReducedKGrid_cP{D}) where D = ntuple(_ -> kG.Ns, D)
gridshape(kG::FullKGrid_cP{D}) where D = ntuple(_ -> kG.Ns, D)

# ---------------------------- BZ to f. irr. BZ -------------------------------

"""
    reduceKGrid(kG::FullKGrid{T}) where T <: cP

Returns the grid on the fully irredrucible BZ.
Filters an arbitrary grid, defined on a full kGrid, so that only 
the lower triangle remains, i.e.
for any (x_1, x_2, ...) the condition x_1 >= x_2 >= x_3 ... 
is fulfilled.
"""
function reduceKGrid(kG::FullKGrid{cP,D}) where D 
    kGrid = reshape(kG.kGrid, gridshape(kG))
    ϵkGrid = reshape(kG.ϵkGrid, gridshape(kG))
    ind = collect(Base.product([1:kG.Ns for Di in 1:D]...))

    ll = floor(Int,size(kGrid,1)/2 + 1) - 1
    la = ceil(Int,kG.Ns/2)
    index = if D == 2
        [[x,y]  for x=la:ll+la for y=la:x]
    elseif D == 3
        [[x,y,z] for x=la:ll+la for y=la:x for z = la:y]
    else
        error("D ∉ [2,3] not implemented yet!")
    end

    ind_red = GridInd{D}(undef, length(index))
    grid_red = GridPoints{D}(undef, length(index))
    ϵk_red = GridDisp(undef, length(index))

    # Compute reduced arrays
    for (i,ti) in enumerate(index)
        ind_red[i] = ind[ti...]
        grid_red[i] = kGrid[ti...]
        ϵk_red[i] = ϵkGrid[ti...]
    end

    kMult, expand_perms = build_expand_mapping_SC(D, kG.Ns, ind_red)
    expand_cache = Array{Complex{Float64}}(undef, gridshape(kG))

    return ReducedKGrid_cP{D}(kG.Nk, kG.Ns, kG.t, ind_red, kMult, grid_red, ϵk_red,
                              expand_perms, expand_cache, kG.fftw_plan)
end

"""
    expandKArr(kG::ReducedKGrid{T1}, arr::Array{T2,1})

Expands array of values on reduced k grid back to full BZ.
"""

function expandKArr(kG::ReducedKGrid_cP{D}, arr::AbstractArray{T, 1})::AbstractArray{T, D} where {T,D}
    length(arr) != length(kG.kInd) && throw(ArgumentError("length of k grid ($(length(kG.kInd))) and argument ($(length(arr))) not matching"))
    res = similar(arr, gridshape(kG)...)
    expandKArr!(kG, res, arr)
    return res
end

function expandKArr!(kG::ReducedKGrid_cP, arr::Array{Complex{Float64}, 1})
    for (ri,perms) in enumerate(kG.expand_perms)
        @simd for p in perms
            @inbounds kG.expand_cache[p] = arr[ri]
        end
    end
end

function expandKArr!(kG::ReducedKGrid_cP, res::AbstractArray{T,D}, arr::Array{T, 1}) where {T,D}
    for (ri,perms) in enumerate(kG.expand_perms)
        @simd for p in perms
            @inbounds res[p] = arr[ri]
        end
    end
end


function reduceKArr(kG::ReducedKGrid_cP{D}, arr::AbstractArray{T,D}) where {T,D}
    res = similar(arr, length(kG.kInd))
    reduceKArr!(kG, res, arr)
    return res
end

function reduceKArr!(kG::ReducedKGrid_cP{D}, res::AbstractArray{T,1}, arr::AbstractArray{T,D}) where {T,D}
    for (i,ki) in enumerate(kG.kInd)
        @inbounds res[i] = arr[ki...]
    end
end

gen_ϵkGrid(::Type{cP}, kGrid::GridPoints, t::T) where T <: Real = collect(map(kᵢ -> -2*t*sum(cos.(kᵢ)), kGrid))
ifft_post(kG::ReducedKGrid_cP, x::Array{T,N}) where {N, T <: Number} = ShiftedArrays.circshift(x, floor.(Int, gridshape(kG) ./ 2) .+ 1)

function build_expand_mapping_SC(D::Int, Ns::Int, ind_red::Array)
    expand_perms = Vector{Vector{CartesianIndex{D}}}(undef, length(ind_red))
    kMult = Array{Int,1}(undef, length(ind_red))

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
