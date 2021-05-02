module Dispersions

using Combinatorics

# general types
export KGrid, FullKGrid, ReducedKGrid

# access functions
export gridPoints, Nk

# grid functions
export reduceKGrid, reduceKArr, reduceKArr_reverse, expandKArr, conv_transform

# grids 
export gen_cP_kGrid, cP_2D, cP_3D

# sum types
export KSum
# functions
export kintegrate


include("Types.jl")
include("SC.jl")
include("BZIntegration.jl")
include("IO_SC.jl")

end
