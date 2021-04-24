# ================================================================================ #
#                                   Type Defs                                      #
# ================================================================================ #

abstract type KGridType end
abstract type KGrid{T <: KGridType} end
abstract type FullKGrid{T} <: KGrid{T} end
abstract type ReducedKGrid{T} <: KGrid{T} end

# --------------------------------- convenience defs -------------------------------
const GridInd2D = Array{Tuple{Int,Int},1}
const GridInd3D = Array{Tuple{Int,Int,Int},1}
const GridPoints2D = Array{Tuple{Float64,Float64},1}
const GridPoints3D = Array{Tuple{Float64,Float64,Float64},1}
const GridDisp = Array{Float64,1}

# Trait based dispatch functionality for GridPoints and GridInd
struct Not{T} end 
struct IsGridPoints{T} end
isGridPoints(::T) where {T}  = Not{IsGridPoint{T}}
isGridPoints(::GridPoints2D) = IsGridPoints{GridPoints2D}
isGridPoints(::GridPoints3D) = IsGridPoints{GridPoints3D}

struct IsGridInd{T} end
isGridInd(::T) where {T}  = Not{IsGridInd{T}}
isGridInd(::GridInd2D) = IsGridInd{GridInd2D}
isGridInd(::GridInd3D) = IsGridInd{GridInd3D}

# ================================================================================ #
#                                   Interface                                      #
# ================================================================================ #

# ------------------------------ Helper Functions -----------------------------

"""
    Nk(kG::T) where T <: KGrid

Total number of k points (length of `kGrid.kGrid` for full grids). 
"""
Nk(kG::T) where T <: KGrid = kG.Nk

"""
    gridPoints(kG::T)::Int where T <: KGrid

Number of grid points as integer value. `Nk` needs to be accessible for all implemented 
full and reduced k grids.
"""
gridPoints(kG::T) where T <: KGrid = kG.kGrid

# ------------------------------ Helper Functions -----------------------------
"""
    reduceKGrid(kGrid::FullKGrid{T}) where T <: KGridType

Returns the grid on the fully irredrucible BZ.
"""
reduceKGrid(kGrid::FullKGrid{T}) where T <: KGridType = throw(MethodError("KGrid Instance not found")) 

expandKGrid(kGrid::FullKGrid{T}, arr::Array) where T <: KGridType = throw(MethodError("KGrid Instance not found")) 


# ------------------------------ Private Functions -----------------------------
"""
    gen_ϵkGrid(::T1, kGrid::T2, t::T3) where {T1 <: KGridType, T2 <: Array, T3 <: Number}

Generates dispersion relation for array of k-points. Usually only accessed by constructor of
KGrid
"""
gen_ϵkGrid(::Type{T1}, kGrid::T2, t::T3) where {T1 <: KGridType, T2 <: Array, T3 <: Number} = throw(MethodError("KGrid Instance not found"))

