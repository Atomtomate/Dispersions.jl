#TODO: implement and test transform_to_first_BZ
#TODO: documentation and header style from LadderDGA.jl 
#
# ================================================================================ #
#                                   Type Defs                                      #
# ================================================================================ #

# --------------------------------- convenience defs -------------------------------
const GridInd{D} = Vector{CartesianIndex{D}}
const GridPoints{D} = Vector{NTuple{D,Float64}}
const GridDisp = Union{Array{Float64,1}, Array{ComplexF64,3}}

# ------------------------------------ Grids -----------------------------------
abstract type KGridType end

# The following functions are expected to be implemented by all grid types.

gen_sampling(gt::GT, D::Int, Nk::Int) where {GT} =
    throw(ArgumentError("Cannot generate sampling! Grid type $gt unkown!"))

basis_transform(gt::GT, v::AbstractVector) where {GT} =
    throw(ArgumentError("Cannot basis transform! Grid type $gt unkown!"))

reduce_KGrid(gt::GT, D::Int, Nk::Int, kGrid::Vector) where {GT} =
    throw(ArgumentError("Cannot reduce k grid! Grid type $gt unkown!"))


# """
#     symmetry_points(::GridType)

# Returns:
# -------------
# PointList : `Vector{GridInd}`, List of high symmetry points, each lattice overloads this.
# """
# function symmetry_points(kG::KGrid)::GridPoints where GT
#     return GridPoints[]
# end

# """
#     symmetry_paths(::GridType, N::Int)

# Returns:
# -------------
# PathsList : `Vector`, List of symmetry paths, each lattice overloads this.
# """
# function symmetry_paths(kG::KGrid)::Vector where GT
#     []
# end
# #gen_ϵkGrid(gt, kGrid, t) =  throw(ArgumentError("Cannot generate dispersion relation! Grid type $gt unkown!"))
