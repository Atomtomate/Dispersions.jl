function gen_kGrid(kg::String, Nk::Int; full=false)
    sp = findfirst("-", kg)[1]
    data = [kg[1:(sp-1)], kg[(sp+1):end]]
    grid = if lowercase(data[1]) == "3dsc"
        FullKGrid_cP(3, Nk, parse(Float64, data[2]))
    elseif lowercase(data[1]) == "2dsc"
        FullKGrid_cP(2, Nk, parse(Float64, data[2]))
    elseif lowercase(data[1]) == "p6m"
        FullKGrid_p6m(Nk, parse(Float64, data[2]))
    elseif lowercase(data[1]) == "file"
        FullKGrid_File(data[2])
    else
        throw(ArgumentError("Unkown grid type"))
    end
    res = full ? grid : reduceKGrid(grid)
    return res
end

# ================================================================================ #
#                                   Interface                                      #
# ================================================================================ #

# ------------------------------ Helper Functions -----------------------------

"""
    gridshape(kG::T) where T <: KGrid

shape of kGrid (e.g. `(kG.Ns, kG.Ns)` for 2D sc) 
"""
gridshape(kG) = throw(ArgumentError("KGrid Instance of $(kG) not found"))

"""
    Nk(kG::T) where T <: KGrid

Total number of k points (length of `kGrid.kGrid` for full grids). 
"""
Nk(kG) = kG.Nk

"""
    gridPoints(kG::T)::Int where T <: KGrid

Number of grid points as integer value. `Nk` needs to be accessible for all implemented 
full and reduced k grids.
"""
gridPoints(kG::T) where T <: KGrid = kG.kGrid

"""
    dispersion(kG::T)::Int where T <: KGrid

Returns dispersion relation of grid.
"""
dispersion(kG::T) where T <: KGrid = kG.ϵkGrid

# ------------------------------ Helper Functions -----------------------------
"""
    reduceKGrid(kGrid::FullKGrid{T}) where T <: KGridType

Returns the grid on the fully irredrucible BZ.
"""
reduceKGrid(kG) = throw(ArgumentError("KGrid Instance of $(kG) not found")) 

expandKArr(kG, arr) = throw(ArgumentError("KGrid Instance of $(kG) not found")) 


# ------------------------------ Convolution Functions -----------------------------

function conv(kG::ReducedKGrid, arr1::AbstractArray{Complex{Float64},1}, arr2::AbstractArray{Complex{Float64},1})
    Nk(kG) == 1 && return arr1 .* arr2
    tmp = reshape(fft(expandKArr(kG, arr1)) .* fft(expandKArr(kG, arr2)), gridshape(kG)) |> ifft 
    ifft_post!(typeof(kG), tmp)
    return reduceKArr(kG, tmp) ./ Nk(kG)
end

function conv!(kG::ReducedKGrid, res::AbstractArray{Complex{Float64},1}, arr1::AbstractArray{Complex{Float64},1}, arr2::AbstractArray{Complex{Float64},1})
    Nk(kG) == 1 && (res[:] = arr1 .* arr2)
    expandKArr!(kG, arr1)
    tmp = fft(kG.expand_cache)
    expandKArr!(kG, arr2)
    fft!(kG.expand_cache)
    kG.expand_cache[:] = kG.expand_cache .* tmp
    AbstractFFTs.ldiv!(kG.expand_cache, kG.fftw_plan, kG.expand_cache)
    ifft_post!(typeof(kG), kG.expand_cache)
    reduceKArr!(kG, res, kG.expand_cache) 
    res[:] = res ./ Nk(kG)
end


function conv_fft1(kG::ReducedKGrid, arr1::AbstractArray{Complex{Float64},1}, arr2::AbstractArray{Complex{Float64}})
    Nk(kG) == 1 && return arr1 .* arr2
    newArr = Array{eltype(arr1),1}(undef, length(kG.kInd))
    conv_fft1!(kG, newArr, arr1, reshape(arr2, gridshape(kG)))
    return newArr
end

function conv_fft1!(kG::ReducedKGrid, res::AbstractArray{Complex{Float64},1}, arr1::AbstractArray{Complex{Float64},1}, arr2::AbstractArray{Complex{Float64}})
    Nk(kG) == 1 && (res[:] = arr1 .* arr2)
    expandKArr!(kG, arr1)
    AbstractFFTs.mul!(kG.expand_cache, kG.fftw_plan, kG.expand_cache)
    @simd for i in 1:length(kG.expand_cache)
        @inbounds kG.expand_cache[i] *= arr2[i] 
    end
    AbstractFFTs.ldiv!(kG.expand_cache, kG.fftw_plan, kG.expand_cache)
    ifft_post!(typeof(kG), kG.expand_cache)
    reduceKArr!(kG, res, kG.expand_cache) 
    @simd for i in 1:length(res)
        @inbounds res[i] /= kG.Nk
    end
end

function conv_fft(kG::ReducedKGrid, arr1::AbstractArray{Complex{Float64}}, arr2::AbstractArray{Complex{Float64}})
    Nk(kG) == 1 && return arr1 .* arr2
    reduceKArr(kG, ifft_post!(typeof(kG), ifft(arr1 .* arr2))) ./ Nk(kG)
end

function conv_fft!(kG::ReducedKGrid, res::AbstractArray{Complex{Float64},1}, arr1::AbstractArray{Complex{Float64}}, arr2::AbstractArray{Complex{Float64}})
    Nk(kG) == 1 && (res[:] = arr1 .* arr2)

    @simd for i in eachindex(kG.expand_cache)
        @inbounds kG.expand_cache[i] = arr1[i] .* arr2[i]
    end
    AbstractFFTs.ldiv!(kG.expand_cache, kG.fftw_plan, kG.expand_cache)
    ifft_post!(typeof(kG), kG.expand_cache)
    reduceKArr!(kG, res, kG.expand_cache) 
    @simd for i in 1:length(res)
        @inbounds res[i] /= kG.Nk
    end
end
