# ================================================================================ #
#                                   Type Defs                                      #
# ================================================================================ #

#TODO: exteremely simplistic stub, needs to be replaced by some actual k space integration later on

abstract type KSpaceIntegrator end

struct KSum <: KSpaceIntegrator end


kintegrate(grid::Nothing, arr::AbstractArray) = arr

function kintegrate(grid::T1, arr::AbstractArray; dim=1, type::T2 = KSum()) where {T1 <: ReducedKGrid, T2 <: KSum}
    size(arr)[dim] != length(grid.kMult) && throw(ArgumentError("Dimension does not seem to be on a k grid! Length is $(size(arr)[dim]) but should be $(length(grid.kMult))."))
    return mapslices(sub_arr -> sum(grid.kMult .* sub_arr)/sum(grid.kMult), arr, dims=dim)
end

function kintegrate(grid::T1, arr::AbstractArray; dim=1, type::T2 = KSum()) where {T1 <: FullKGrid, T2 <: KSum}
    !all(size(grid.ϵkGrid) .== size(arr)) && throw(ArgumentError("Array does not seem to be on a k grid! Size is $(size(arr)) but should be $(size(grid.ϵkGrid))."))
    return sum(arr)/length(arr)
end
