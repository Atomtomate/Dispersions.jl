# ================================================================================ #
#                                   Type Defs                                      #
# ================================================================================ #

#TODO: exteremely simplistic stub, needs to be replaced by some actual k space integration later on

abstract type KSpaceIntegrator end

struct KSum <: KSpaceIntegrator end


kintegrate(kG::Nothing, arr::AbstractArray) = arr

function kintegrate(kG::T1, arr::AbstractArray, dim::Int, type::T2 = KSum()) where {T1 <: ReducedKGrid, T2 <: KSum}
    size(arr)[dim] != length(kG.kMult) && throw(ArgumentError("Dimension does not seem to be on a k grid! Length is $(size(arr)[dim]) but should be $(length(kG.kMult))."))
    return mapslices(sub_arr -> kintegrate(kG, sub_arr, type=type), arr, dims=dim)
end

function kintegrate(kG::T1, arr::AbstractArray{T2,1}; type::T3 = KSum()) where {T1 <: ReducedKGrid, T2 <: Number, T3 <: KSum}
    return dot(kG.kMult, arr)/Nk(kG)
end

#function kintegrate(grid::T1, arr::AbstractArray; dim=1, type::T2 = KSum()) where {T1 <: FullKGrid, T2 <: KSum}
#    !all(size(grid.ϵkGrid) .== size(arr)) && throw(ArgumentError("Array does not seem to be on a k grid! Size is $(size(arr)) but should be $(size(grid.ϵkGrid))."))
#    return mapslices(sub_arr -> sum(grid.kMult .* sub_arr)/sum(grid.kMult), arr, dims=dim)
#end
