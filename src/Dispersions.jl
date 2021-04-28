module Dispersions

using Combinatorics

# general types
export KGrid, FullKGrid, ReducedKGrid

# access functions
export gridPoints, Nk

# grid functions
export reduceKGrid, reduceKArr, expandKArr

# grids 
export gen_cP_kGrid

# sum types
export KSum
# functions
export kintegrate


include("Types.jl")
include("SC.jl")
include("BZIntegration.jl")
include("IO_SC.jl")

end
