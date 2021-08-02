module Dispersions

using Combinatorics
using AbstractFFTs, FFTW
using HDF5

# general types
export KGrid, FullKGrid, ReducedKGrid

# access functions
export gridPoints, Nk, gridshape, dispersion

# grid functions
export reduceKArr, reduceKArr!, expandKArr, expandKArr!, conv, conv!, conv_fft, conv_fft!, conv_fft1, conv_fft1!

# grids 
export gen_kGrid, cP_2D, cP_3D, p6m, FileDisp

# sum types
export KSum
# functions
export kintegrate


include("Types.jl")
include("common.jl")
include("SC.jl")
include("hexagonal.jl")
include("fileDisp.jl")
include("BZIntegration.jl")
include("IO_SC.jl")

end

#TODO:iplement input fro mfile
#TODO: implement https://arxiv.org/pdf/2104.05856.pdf
