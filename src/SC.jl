import Base.collect


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
        return FullKGrid_cP_2D(Nk, t)
    elseif D == 3
        return FullKGrid_cP_3D(Nk, t)
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

abstract type cP_2D <: KGridType end
"""
    FullKGrid_cP_2D  <: FullKGrid{cP_2D}

Fields
-------------
- **`kGrid`** : `Array{Tuple{Float64, ...}}` of kGrids. Each element is a D-tuple
"""
struct FullKGrid_cP_2D  <: FullKGrid{cP_2D}
    Nk::Int
    Ns::Int
    kGrid::GridPoints2D
    ϵkGrid::GridDisp
    t::Float64
    function FullKGrid_cP_2D(Nk::Int, t::Float64)
        #isodd(Nk) && throw(ArgumentError("Only implemented for even k grids"))
        kx = [(2*π/Nk) * j - π for j in 1:Nk]
        kGrid  = collect(Base.product([kx for Di in 1:2]...))[:]
        new(Nk^2, Nk, kGrid, gen_ϵkGrid(cP_2D,kGrid,t),t)
    end
end

"""
    ReducedKGrid_cP_2D  <: ReducedKGrid{cP_2D}

Reduced k grid only containing points necessary for computation of quantities involving 
this grid type and multiplicity of each stored point.

Fields
-------------
- **`Nk`**    : `Int` k points in full grid.
- **`kInd`**  : `Array{Tuple{Int,...}}` indices in full grid. used for reconstruction of full grid.
- **`kMult`** : `Array{Float64,1}` multiplicity of point, used for calculations involving reduced k grids.
- **`kGrid`** : `Array{Tuple{Float64,...}}` k points of reduced grid.
"""
struct ReducedKGrid_cP_2D  <: ReducedKGrid{cP_2D}
    Nk::Int
    kInd::GridInd2D
    kMult::Array{Float64,1}
    kGrid::GridPoints2D
    ϵkGrid::GridDisp
end

# -------------------------------------------------------------------------------- #
#                                   Interface                                      #
# -------------------------------------------------------------------------------- #

# ---------------------------- BZ to f. irr. BZ -------------------------------
gen_ϵkGrid(::Type{cP_2D}, kGrid::GridPoints2D, t::T1) where T1 <: Number = collect(map(kᵢ -> t*(cos(kᵢ[1])+cos(kᵢ[2])), kGrid))

"""
    reduceKGrid(kG::FullKGrid{T}) where T <: Union{cP_2D, cP_3D}

Returns the grid on the fully irredrucible BZ.
Filters an arbitrary grid, defined on a full kGrid, so that only 
the lower triangle remains, i.e.
for any (x_1, x_2, ...) the condition x_1 >= x_2 >= x_3 ... 
is fulfilled.
"""
function reduceKGrid(kG::FullKGrid{cP_2D})
    s = [kG.Ns for i in 1:2]
    grid_tmp = cut_mirror(reshape(kG.kGrid, s...))
    ϵk_tmp = cut_mirror(reshape(kG.ϵkGrid, s...))
    ind = cut_mirror(collect(Base.product([1:kG.Ns for Di in 1:2]...)))

    index = [[x,y] for x=1:size(ind,1) for y=1:x]

    ind_red = GridInd2D(undef, length(index))
    grid_red = GridPoints2D(undef, length(index))
    ϵk_red = GridDisp(undef, length(index))

    for (i,ti) in enumerate(index)
        ind_red[i] = ind[ti...]
        grid_red[i] = grid_tmp[ti...]
        ϵk_red[i] = ϵk_tmp[ti...]
    end
	kmult = kGrid_multiplicity_cP(ind_red)
    return ReducedKGrid_cP_2D(kG.Nk, ind_red, kmult, grid_red, ϵk_red)
end


"""
    expandKArr(kG::ReducedKGrid{T1}, arr::Array{T2,1}) where {T1 <: Union{cP_2D,cP_3D}, T2 <: Any

Expands array of values on reduced k grid back to full BZ.
"""
function expandKArr(kG::ReducedKGrid{cP_2D}, arr::Array{T, 1}) where T
    length(arr) != length(kG.kInd) && throw(ArgumentError("length of k grid ($(length(kG.kInd))) and argument ($(length(arr))) not matching"))
    N = maximum(maximum.(kG.kInd))
    newArr = Array{eltype(arr)}(undef, (N*ones(Int64, 2))...)
    for (ri, redInd) in enumerate(kG.kInd)
        perms = unique(collect(permutations(redInd)))
        for p in perms
            newArr[p...] = arr[ri]
        end
    end
    minimum(minimum.(kG.kInd)) > 1 && expand_mirror!(newArr)
    return newArr
end

"""
    kGrid_multiplicity(kIndices)

Compute multiplicity for each k point over a given index
array of a reduced kGrid.
"""
#function multiplicity(kIndices::ReducedKIndices_cP)
#    min_ind = minimum(kIndices.grid)
#    max_ind = maximum(kIndices.grid)
#    borderFactor(el) = prod((el[i] == min_ind[i] || el[i] == max_ind[i]) ? 0.5 : 1.0 for i in 1:length(min_ind))
#    length(min_ind) == 2 ? map(el -> borderFactor(el)*8/((el[2]==el[1]) + 1), kIndices.grid) :
#            map(el -> borderFactor(el)*48/( (el[2]==el[1]) + (el[3]==el[2]) + 3*(el[3]==el[1]) + 1), kIndices.grid)
#end
# ================================================================================ #
#                                Simple Cubic 3D                                   #
# ================================================================================ #

# -------------------------------------------------------------------------------- #
#                                     Types                                        #
# -------------------------------------------------------------------------------- #
abstract type cP_3D <: KGridType end

"""
    FullKGrid_cP_3D  <: FullKGrid{cP_3D}

Fields
-------------
- **`kGrid`** : `Array{Tuple{Float64, ...}}` of kGrids. Each element is a D-tuple
"""
struct FullKGrid_cP_3D  <: FullKGrid{cP_3D}
    Nk::Int
    Ns::Int
    kGrid::GridPoints3D
    ϵkGrid::GridDisp
    t::Float64
    function FullKGrid_cP_3D(Nk::Int, t::Float64)
        kx = [(2*π/Nk) * j - π for j in 1:Nk];
        kGrid  = collect(Base.product([kx for Di in 1:3]...))[:]
        new(Nk^3, Nk, kGrid, gen_ϵkGrid(cP_3D,kGrid,t),t)
    end
end


"""
    ReducedKGrid_cP_3D  <: ReducedKGrid{cP_3D}

Reduced k grid only containing points necessary for computation of quantities involving 
this grid type and multiplicity of each stored point.

Fields
-------------
- **`Nk`**    : `Int` k points in full grid.
- **`kInd`**  : `Array{Tuple{Int,...}}` indices in full grid. used for reconstruction of full grid.
- **`kMult`** : `Array{Float64,1}` multiplicity of point, used for calculations involving reduced k grids.
- **`kGrid`** : `Array{Tuple{Float64,...}}` k points of reduced grid.
"""
struct ReducedKGrid_cP_3D  <: ReducedKGrid{cP_3D}
    Nk::Int
    kInd::GridInd3D
    kMult::Array{Float64,1}
    kGrid::GridPoints3D
    ϵkGrid::GridDisp
end


function reduceKGrid(kG::FullKGrid{cP_3D})
    s = [kG.Ns for i in 1:3]
    grid_tmp = cut_mirror(reshape(kG.kGrid, s...))
    ϵk_tmp = cut_mirror(reshape(kG.ϵkGrid, s...))
    ind = cut_mirror(collect(Base.product([1:kG.Ns for Di in 1:3]...)))

    index = [[x,y,z] for x=1:size(ind,1) for y=1:x for z = 1:y]
    ind_red = GridInd3D(undef, length(index))
    grid_red = GridPoints3D(undef, length(index))
    ϵk_red = GridDisp(undef, length(index))

    for (i,ti) in enumerate(index)
        ind_red[i] = ind[ti...]
        grid_red[i] = grid_tmp[ti...]
        ϵk_red[i] = ϵk_tmp[ti...]
    end
	kmult = kGrid_multiplicity_cP(ind_red)
    return ReducedKGrid_cP_3D(kG.Nk, ind_red, kmult, grid_red, ϵk_red)
end


function expandKArr(kG::ReducedKGrid{cP_3D}, arr::Array)
    length(arr) != length(kG.kInd) && throw(ArgumentError("length of k grid ($(length(kG.kInd))) and argument ($(length(arr))) not matching"))
    N = maximum(maximum.(kG.kInd))
    newArr = Array{eltype(arr)}(undef, (N*ones(Int64, 3))...)
    for (ri, redInd) in enumerate(kG.kInd)
        perms = unique(collect(permutations(redInd)))
        for p in perms
            newArr[p...] = arr[ri]
        end
    end
    expand_mirror!(newArr)
    return newArr
end


# ---------------------------- Auxilliary Functions -------------------------------
# TODO: avoid reshape
@inbounds cut_mirror(arr::Base.Iterators.ProductIterator) = cut_mirror(collect(arr))
@inbounds cut_mirror(arr::Array{T, 2}) where T = arr[ceil(Int,size(arr,1)/2):end, ceil(Int,size(arr,2)/2):end]
@inbounds cut_mirror(arr::Array{T, 3}) where T = arr[ceil(Int,size(arr,1)/2):end, ceil(Int,size(arr,2)/2):end, ceil(Int,size(arr,3)/2):end]
#reverse cut. This is a helper function to avoid reversing the array after fft-convolution trick. assumes reserved input and returns correct array including cut
@inbounds ifft_cut_mirror(arr::Base.Iterators.ProductIterator) = ifft_cut_mirror(collect(arr))
@inbounds ifft_cut_mirror(arr::Array{T, 2}) where T = arr[end-1:-1:Int64(size(arr,1)/2-1), 
                                                          end-1:-1:Int64(size(arr,2)/2-1)]
@inbounds ifft_cut_mirror(arr::Array{T, 3}) where T = arr[end-1:-1:Int64(size(arr,1)/2-1), 
                                                          end-1:-1:Int64(size(arr,2)/2-1), 
                                                          end-1:-1:Int64(size(arr,3)/2-1)]


function expand_mirror!(arr::Array{T, 2}) where T 
    N = size(arr,1)
    al = ceil(Int,N/2)

    arr[1:al-1,al:end] = arr[end-1:-1:iseven(N)+al,al:end]
    arr[1:end,1:al-1] = arr[1:end,end-1:-1:al+iseven(N)]

    return arr
end

function expand_mirror!(arr::Array{T, 3}) where T
    N = size(arr,1)
    al = ceil(Int,N/2)

    Nl = N - al
    # copy values along anti-diagonals
    for ll in al:N 
        for i in 0:Nl-1
            for j in 1:N-1
                arr[ll,mod1(al-j,N),mod1(al+j+i,N)] = arr[ll,mod1(al,N),mod1(al+i,N)]
            end
        end
        for i in 0:Nl
            for j in 1:(N-i-1)
                arr[ll,N-i-j,j] = arr[ll,N-i,N]
            end
        end
    end
    for ll in al-1:-1:1
        arr[ll,:,:] = circshift(arr[ll+1,:,:],(1,0))
    end
    return arr
end


"""
	kGrid_multiplicity_cP(kIndices)

Given a set of reduced indices, produce list of multiplicities for each point
"""
function kGrid_multiplicity_cP(kIndices)
    min_ind = minimum(kIndices)
    max_ind = maximum(kIndices)
    function borderFactor(el)
        val = 1.0
        for i in 1:length(el)
            val = if (el[i] == min_ind[i] || el[i] == max_ind[i]) val*0.5 else val end
        end
        return val
    end
    if length(min_ind) == 2
        res = map(el -> borderFactor(el)*8/((el[2]==el[1]) + 1), kIndices)
    elseif length(min_ind) == 3
        res = map(el -> borderFactor(el)*48/( (el[2]==el[1]) + (el[3]==el[2]) + 3*(el[3]==el[1]) + 1), kIndices)
    else
        throw("Multiplicity of k points only implemented for 2D and 3D")
    end
    return res
end

gen_ϵkGrid(::Type{cP_3D}, kGrid::GridPoints3D, t::T1) where T1 <: Number = collect(map(kᵢ -> t*(cos(kᵢ[1]+kᵢ[2]+kᵢ[3])), kGrid))
