# ================================================================================ #
#                                   Type Defs                                      #
# ================================================================================ #

# ------------------------------------ Grids -----------------------------------
abstract type KGridType end
abstract type KGrid{T <: KGridType, D} <: KGridType end
abstract type FullKGrid{T <: KGridType, D} <: KGrid{T, D} end
abstract type ReducedKGrid{T <: KGridType, D} <: KGrid{T, D} end

# --------------------------------- convenience defs -------------------------------
const GridInd{D} = Array{NTuple{D,Int}}
const GridPoints{D} = Array{NTuple{D,Float64}}
const GridDisp = Array{Float64,1}
