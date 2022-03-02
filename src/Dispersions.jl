module Dispersions

using Combinatorics
using AbstractFFTs, FFTW
using LinearAlgebra
using ShiftedArrays
using FastGaussQuadrature

# general types
export KGrid

# access functions
export gridPoints, Nk, gridshape, dispersion

# grid functions
export conv, conv!, conv_fft, conv_fft!, conv_fft1, conv_fft1!

# grids 
export gen_kGrid, SC, FCC

# sum types
export KSum
# functions
export kintegrate


include("Types.jl")
include("KGrid.jl")
include("common.jl")
include("IO.jl")
include("BZIntegration.jl")


end

#TODO: implement input from basis vector
#TODO: implement https://arxiv.org/pdf/2104.05856.pdf
