# ================================================================================ #
#                                   Type Defs                                      #
# ================================================================================ #

#TODO: exteremely simplistic stub, needs to be replaced by some actual k space integration later on

abstract type KSpaceIntegrator end

struct KSum <: KSpaceIntegrator end


kintegrate(grid::Nothing, arr::AbstractArray) = arr

function kintegrate(grid::T1, arr::AbstractArray; type::T2 = KSum()) where {T1 <: ReducedKGrid, T2 <: KSum}
    return sum(grid.kMult .* arr)/sum(grid.kMult)
end

function kintegrate(grid::T1, arr::AbstractArray; type::T2 = KSum()) where {T1 <: FullKGrid, T2 <: KSum}
    return sum(arr)/length(arr)
end
