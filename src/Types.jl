# ================================================================================ #
#                                   Type Defs                                      #
# ================================================================================ #

# ------------------------------------ Grids -----------------------------------
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


