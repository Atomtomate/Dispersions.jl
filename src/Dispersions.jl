module Dispersions

using Combinatorics
using AbstractFFTs, FFTW
using LinearAlgebra
using ShiftedArrays
using FastGaussQuadrature
using EquivalenceClassesConstructor

# general types
export KGrid

# access functions
export gridPoints, Nk, gridshape, dispersion

# grid functions
export reduceKArr,reduceKArr!,expandKArr,expandKArr!,
    conv,conv!,conv_fft,conv_fft!,conv_fft1,conv_fft1!
export conv_noPlan,conv_noPlan!,conv_fft_noPlan,conv_fft_noPlan!,conv_fft1_noPlan,conv_fft1_noPlan!
    

# grids 
export gen_kGrid, cP, cF, p6m

# sum types
export KSum
# functions
export kintegrate


include("Types.jl")
include("KGrid.jl")
include("common.jl")
include("cP.jl")
include("cF.jl")
include("hexagonal.jl")
include("BZIntegration.jl")
include("IO.jl")


end

#TODO: implement input from file
#TODO: implement https://arxiv.org/pdf/2104.05856.pdf
