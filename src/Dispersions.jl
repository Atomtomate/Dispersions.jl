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
export gen_sampling, transform_to_first_BZ, gridPoints, Nk, gridshape, dispersion, grid_type, grid_dimension, grid_dimension, ϵ_k_plus_q

# grid functions
export reduceKArr,reduceKArr!,expandKArr,expandKArr!,
    conv,conv!,conv_fft,conv_fft!,conv_fft1,conv_fft1!
export conv_noPlan,conv_noPlan!,conv_fft_noPlan,conv_fft_noPlan!,conv_fft1_noPlan,conv_fft1_noPlan!
    
# helper functions
export build_q_lookup

# grids 
export gen_kGrid, cP, cPnn, cF, cI, Hofstadter

# sum types
export KSum

# integrate functions
export kintegrate


include("Types.jl")
include("KGrid.jl")
include("common.jl")
include("cP.jl")
include("cPnn.jl")
include("cPnn_mag.jl")
include("cF.jl")
include("cI.jl")
include("BZIntegration.jl")
include("IO.jl")

include("helpers.jl")


end

#TODO: implement input from file
#TODO: implement https://arxiv.org/pdf/2104.05856.pdf
