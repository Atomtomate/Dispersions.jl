import Base.collect

# ================================================================================ #
#                                    Exports                                       #
# ================================================================================ #

export gen_cP_kGrid, reduce, reduce_kGrid_cP, expand


# ================================================================================ #
#                               Type Definitions                                   #
# ================================================================================ #
abstract type KIndices{T} end
abstract type FullKIndices{T <: Base.Iterators.ProductIterator} <: KIndices{T} end
abstract type ReducedKIndices{T <: Array} <: KIndices{T} end

abstract type KPoints{T} end
abstract type FullKPoints{T <: Base.Iterators.ProductIterator} <: KPoints{T} end
abstract type ReducedKPoints{T} <: KPoints{T} end

abstract type KGrid{T1 <: KIndices, T2 <: KPoints} end

abstract type Dispersion{T <: KGrid} end

# -------------------------------------------------------------------------------- #
#                                 Simple Cubic                                     #
# -------------------------------------------------------------------------------- #

# ------------------------------- Type Defs --------------------------------
abstract type KGrid_cP{T1, T2} <: KGrid{T1, T2} end


# ------------------------------- Full Grids -------------------------------
const Ind2D = Tuple{UnitRange{Int64},UnitRange{Int64}}
const Ind3D = Tuple{UnitRange{Int64},UnitRange{Int64},UnitRange{Int64}}
const GridP2D = Tuple{Array{Float64,1},Array{Float64,1}}
const GridP3D = Tuple{Array{Float64,1},Array{Float64,1},Array{Float64,1}}

struct FullKIndices_cP{Ind<:Union{Ind2D, Ind3D}} <: FullKIndices{Base.Iterators.ProductIterator{Ind}} 
    grid::Base.Iterators.ProductIterator{Ind}
end
collect(indices::FullKIndices_cP) = collect(indices.grid)

struct FullKPoints_cP{gInd<:Union{GridP2D, GridP3D}} <: FullKPoints{Base.Iterators.ProductIterator{gInd}} 
    grid::Base.Iterators.ProductIterator{gInd}
end
collect(grid::FullKPoints_cP) = collect(grid.grid)

struct FullKGrid_cP{IndT <: FullKIndices_cP, GridT <: FullKPoints_cP}  <: KGrid_cP{IndT, GridT}
    Nk::Int64
    kInd::IndT
    kGrid::GridT
    FullKGrid_cP{IndT, GridT}(Nk::Int64) where {IndT <: FullKIndices_cP, GridT <: FullKPoints_cP} = (
                               kx = [(2*π/Nk) * j - π for j in 1:Nk];
                               D = IndT == FullKIndices_cP{Ind2D} ? 2 : 3;
                               ind = IndT(Base.product([1:(Nk) for Di in 1:D]...));
                               kGrid  = GridT(Base.product([kx for Di in 1:D]...));
                               new{IndT, GridT}(Nk, ind, kGrid)
                              )
end

# ------------------------------- Reduced Grids -------------------------------
const rInd2D = Tuple{Int64,Int64}
const rInd3D = Tuple{Int64,Int64,Int64}
const rGridP2D = Tuple{Float64,Float64}
const rGridP3D = Tuple{Float64,Float64,Float64}

struct ReducedKIndices_cP{Ind<:Union{rInd2D, rInd3D}} <: ReducedKIndices{Array{Ind,1}} 
    grid::Array{Ind,1}
end

struct ReducedKPoints_cP{gInd<:Union{rGridP2D, rGridP3D}} <: ReducedKPoints{Array{gInd,1}} 
    grid::Array{gInd,1}
end


struct ReducedKGrid_cP{IndT <: ReducedKIndices_cP, GridT <: ReducedKPoints_cP}  <: KGrid_cP{IndT, GridT}
    Nk::Int64
    kInd::IndT
    kMult::Array{Float64,1}
    kGrid::GridT
end

# -------------------------------- Grid Data  ---------------------------------
abstract type GridData{DataType <: Number, GridType <: Any} end

#= struct GridData_cP <: GridData =#
#=     grid::Ref{} =#
#= end =#

#= struct GridData_red_cP <: GridData =#
#=     grid::Ref{} =#
#= end =#

# ================================================================================ #
#                                   Functions                                      #
# ================================================================================ #

# ------------------------------ Helper Functions -----------------------------
indices(grid::KGrid_cP) = grid.kInd
gridPoints(grid::KGrid_cP) = grid.kGrid
Nk(grid::KGrid_cP) = grid.Nk


# ---------------------------- BZ to f. irr. BZ -------------------------------

"""
    gen_cP_kGrid(Nk::Int64, D::Int64)

Generates a simple cubic lattice in `D` Dimensions
This can be collected to reduce into a `Nk` times `Nk` [times `Nk`] array, containing
tuples of length `D`.

# Examples
```
julia> gen_cP_kGrid(2, 2)
4,1}}}}(2, Dispersions.FullKIndices_cP{Tuple{UnitRange{Int64},UnitRange{Int64}}}(Base.Iterators.ProductIterator{Tuple{UnitRange{Int64},UnitRange{Int64}}}((1:2, 1:2))), Dispersions.FullKPoints_cP{Tuple{Array{Float64,1},Array{Float64,1}}}(Base.Iterators.ProductIterator{Tuple{Array{Float64,1},Array{Float64,1}}}(([0.0, 3.141592653589793], [0.0, 3.141592653589793]))))
```
"""
function gen_cP_kGrid(Nk::Int64, D::Int64)
    if D == 2
        return FullKGrid_cP{FullKIndices_cP{Ind2D},FullKPoints_cP{GridP2D}}(Nk)
    elseif D == 3
        return FullKGrid_cP{FullKIndices_cP{Ind3D},FullKPoints_cP{GridP3D}}(Nk)
    else
        throw("Simple Cubic ponly implemented for 2D and 3D")
    end
end



"""
    reduce(kGrid::FullKGrid_cP) 

Wrapper function to reduce FullkGrid to reduced kGrid.
"""
function reduce(kGrid::FullKGrid_cP) 
    rInd = kGrid |> indices |> collect |> reduce_kGrid_cP |> ReducedKIndices_cP
    rGrid = kGrid |> gridPoints |> collect |> reduce_kGrid_cP |> ReducedKPoints_cP
    ReducedKGrid_cP(Nk(kGrid), rInd, multiplicity(rInd), rGrid)
end

"""
    reduce_kGrid_cP(kGrid) 

Returns the grid on the fully irredrucible BZ.
Filters an arbitrary grid, defined on a full kGrid, so that only 
the lower triangle remains, i.e.
for any (x_1, x_2, ...) the condition x_1 >= x_2 >= x_3 ... 
is fulfilled.
"""
function reduce_kGrid_cP(kGrid::Union{Array{T,2},Array{T,3}}) where T
    Nk = size(kGrid,1)
    index = ndims(kGrid) == 2 ? [[x,y] for x=1:Nk for y=1:x] : [[x,y,z] for x=1:Nk for y=1:x for z = 1:y]
    grid_red = Array{eltype(kGrid)}(undef, length(index))
    for (i,ti) in enumerate(index)
        grid_red[i] = kGrid[ti...]
    end
    return grid_red
end


"""
    expand_kGrid(reducedInd, reducedArr)

Expands arbitrary Array on reduced k-Grid back to full grid.
This includes restoration of mirror symmetry if the index array 
reducedInd indicates uncomplete grid (by having no (1,1,1) index).
#Examples
```
    mapslices(x->expand_kGrid(qIndices, x, simParams.Nk),sdata(bubble), dims=[2])
```
"""
function expand_kGrid(kGrid, reducedArr::Array{T,1}) where T
#reducedInd::Union{Array{rInd3D,1},Array{rInd2D,1}},
    D = length(reducedInd[1])
    Nk = maximum(maximum.(reducedInd))
    newArr = Array{eltype(reducedArr)}(undef, (Nk*ones(Int64, D))...)
    for (ri,redInd) in enumerate(reducedInd)
        perms = unique(collect(permutations(redInd)))
        for p in perms
            newArr[p...] = reducedArr[ri]
        end
    end
    minimum(minimum.(reducedInd)) > 1 && expand_mirror!(newArr)
    return newArr
end

"""
    kGrid_multiplicity(kIndices)

Compute multiplicity for each k point over a given index
array of a reduced kGrid.
"""
function multiplicity(kIndices::ReducedKIndices_cP)
    min_ind = minimum(kIndices.grid)
    max_ind = maximum(kIndices.grid)
    borderFactor(el) = prod((el[i] == min_ind[i] || el[i] == max_ind[i]) ? 0.5 : 1.0 for i in 1:length(min_ind))
    length(min_ind) == 2 ? map(el -> borderFactor(el)*8/((el[2]==el[1]) + 1), kIndices.grid) :
            map(el -> borderFactor(el)*48/( (el[2]==el[1]) + (el[3]==el[2]) + 3*(el[3]==el[1]) + 1), kIndices.grid)
end

# ---------------------------- Mirror Symmetry -------------------------------

@inbounds cut_mirror(arr::Base.Iterators.ProductIterator) = cut_mirror(collect(arr))
@inbounds cut_mirror(arr::Array{T, 2}) where T = arr[Int(size(arr,1)/2):end, Int(size(arr,2)/2):end]
@inbounds cut_mirror(arr::Array{T, 3}) where T = arr[Int(size(arr,1)/2):end, Int(size(arr,2)/2):end, Int(size(arr,3)/2):end]
#reverse cut. This is a helper function to avoid reversing the array after fft-convolution trick. assumes reserved input and returns correct array including cut
@inbounds ifft_cut_mirror(arr::Base.Iterators.ProductIterator) = ifft_cut_mirror(collect(arr))
@inbounds ifft_cut_mirror(arr::Array{T, 2}) where T = arr[end-1:-1:Int64(size(arr,1)/2-1), 
                                                          end-1:-1:Int64(size(arr,2)/2-1)]
@inbounds ifft_cut_mirror(arr::Array{T, 3}) where T = arr[end-1:-1:Int64(size(arr,1)/2-1), 
                                                          end-1:-1:Int64(size(arr,2)/2-1), 
                                                          end-1:-1:Int64(size(arr,3)/2-1)]


function expand_mirror!(arr::Array{T, 2}) where T <: Any
    al = Int(size(arr,1)/2) - 1

    arr[1:al,al+1:end] = arr[end-1:-1:al+2,al+1:end]
    arr[1:end,1:al] = arr[1:end,end-1:-1:al+2,]
    for i in 1:al
        arr[i,i] = arr[end-i,end-i]
    end
end

function expand_mirror(arr::Array{T, 2}) where T <: Any
    al = size(arr,1) - 2
    size_new = size(arr) .+ al
    res = Array{T}(undef, size_new...)

    res[al+1:end, al+1:end] = arr
    expand_mirror!(res)
    return res
end

function expand_mirror!(arr::Array{T, 3}) where T <: Any
    al = Int(size(arr,1)/2) - 1

    arr[1:al,al+1:end,al+1:end] = arr[end-1:-1:al+2,al+1:end,al+1:end]
    arr[1:end,1:al,al+1:end] = arr[1:end,end-1:-1:al+2,al+1:end]
    for i in 1:al
        arr[i,i,al+1:end] = arr[end-i,end-i,al+1:end]
    end
    arr[1:end,1:end,1:al] .= arr[1:end,1:end,end-1:-1:al+2]
    for i in 1:al
        arr[i,i,i] = arr[end-i,end-i,end-i]
    end
end

function expand_mirror(arr::Array{T, 3}) where T <: Any
    add_length = size(arr,1) - 2
    size_new = size(arr) .+ add_length
    res = Array{T}(undef, size_new...)

    res[add_length:end, add_length:end, add_length:end] = arr
    expand_mirror!(res)
    return res
end
# ================================================================================ #
#                                  Dispersions                                     #
# ================================================================================ #


"""
    squareLattice_ek_grid(kgrid)

Computes 0.5 [cos(k_x) + ... + cos(k_D)] and returns a grid with Nk points.
"""
squareLattice_ekGrid(kgrid)  = ((length(first(kgrid)) == 3 ? -0.40824829046386301636 : -0.5) * 
                                 sum([cos(kᵢ) for kᵢ in k]) for k in kgrid)


function gen_squareLattice_ekq_grid(kList::Any, qList::Any)
    gen_squareLattice_ekq_grid(collect.(kList), collect.(qList))
end

function gen_squareLattice_ekq_grid(kList::Array, qList::Array)
    tsc =  length(first(kList)) == 3 ? -0.40824829046386301636 : -0.5
    res = zeros(length(kList),length(qList))
    for (ki,k) in enumerate(kList)
        for (qi,q) in enumerate(qList)
            @inbounds res[ki,qi] = tsc.*sum(cos.(k .+ q))
        end
    end
    return res
end

#TODO: generalize to 3D, better abstraction
function gen_squareLattice_full_ekq_grid(kList::Array{Tuple{Float64,Float64},1}, qList::Array{Tuple{Float64,Float64},1})
    res = zeros(length(kList),length(qList), 8) # There are 8 additional 
    tsc = -0.5
    for (ki,k) in enumerate(kList)
        for (qi,q) in enumerate(qList)
            @inbounds res[ki,qi,1] = tsc*(cos(k[1] + q[1]) + cos(k[2] + q[2]))
            @inbounds res[ki,qi,2] = tsc*(cos(k[1] + q[1]) + cos(k[2] - q[2]))
            @inbounds res[ki,qi,3] = tsc*(cos(k[1] - q[1]) + cos(k[2] + q[2]))
            @inbounds res[ki,qi,4] = tsc*(cos(k[1] - q[1]) + cos(k[2] - q[2]))
            @inbounds res[ki,qi,5] = tsc*(cos(k[1] + q[2]) + cos(k[2] + q[1]))
            @inbounds res[ki,qi,6] = tsc*(cos(k[1] + q[2]) + cos(k[2] - q[1]))
            @inbounds res[ki,qi,7] = tsc*(cos(k[1] - q[2]) + cos(k[2] + q[1]))
            @inbounds res[ki,qi,8] = tsc*(cos(k[1] - q[2]) + cos(k[2] - q[1]))
        end
    end
    return res
end

function gen_squareLattice_full_ekq_grid(kList::Array{Tuple{Float64,Float64,Float64},1}, qList::Array{Tuple{Float64,Float64,Float64},1})
    perm = permutations([1,2,3])
    res = zeros(length(kList),length(qList), 8*length(perm)) # There are 8 additional 
    tsc = -0.40824829046386301636
    for (ki,k) in enumerate(kList)
        for (qi,q) in enumerate(qList)
            ind = 0
            for p in perm
                res[ki,qi,ind+1] = tsc*(cos(k[p[1]] + q[1]) + cos(k[p[2]] + q[2]) + cos(k[p[3]] + q[3]))
                res[ki,qi,ind+2] = tsc*(cos(k[p[1]] + q[1]) + cos(k[p[2]] + q[2]) + cos(k[p[3]] - q[3]))
                res[ki,qi,ind+3] = tsc*(cos(k[p[1]] + q[1]) + cos(k[p[2]] - q[2]) + cos(k[p[3]] + q[3]))
                res[ki,qi,ind+4] = tsc*(cos(k[p[1]] + q[1]) + cos(k[p[2]] - q[2]) + cos(k[p[3]] - q[3]))
                res[ki,qi,ind+5] = tsc*(cos(k[p[1]] - q[1]) + cos(k[p[2]] + q[2]) + cos(k[p[3]] + q[3]))
                res[ki,qi,ind+6] = tsc*(cos(k[p[1]] - q[1]) + cos(k[p[2]] + q[2]) + cos(k[p[3]] - q[3]))
                res[ki,qi,ind+7] = tsc*(cos(k[p[1]] - q[1]) + cos(k[p[2]] - q[2]) + cos(k[p[3]] + q[3]))
                res[ki,qi,ind+8] = tsc*(cos(k[p[1]] - q[1]) + cos(k[p[2]] - q[2]) + cos(k[p[3]] - q[3]))
                ind += 8
            end
        end
    end
    return res
end

