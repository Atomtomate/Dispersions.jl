function gen_kGrid(kg::String, Nk::Int)
    sp = findfirst("-", kg)[1]
    data = [kg[1:(sp-1)], kg[(sp+1):end]]
    grid = if lowercase(data[1]) == "3dsc"
        KGrid(SC, 3, Nk, parse(Float64, data[2]))
    elseif lowercase(data[1]) == "2dsc"
        KGrid(SC, 2, Nk, parse(Float64, data[2]))
    elseif lowercase(data[1]) == "fcc"
        KGrid_cF(Nk, parse(Float64, data[2]))
    elseif lowercase(data[1]) == "p6m"
        KGrid_p6m(Nk, parse(Float64, data[2]))
    else
        throw(ArgumentError("Unkown grid type"))
    end
    return grid
end



# ================================================================================ #
#                                   Interface                                      #
# ================================================================================ #
export conv, conv!, conv_fft, conv_fft!, conv_fft1, conv_fft1!

#TODO: assumed fields: kInd, expand_erms, expand_cache, ϵkGrid, Nk, kGrid

# ------------------------------ Helper Functions -----------------------------

"""
    gridshape(kG::T) where T <: KGrid

shape of kGrid (e.g. `(kG.Ns, kG.Ns)` for 2D sc) 
"""
gridshape(kG) = throw(ArgumentError("KGrid Instance of $(typeof(kG)) not found"))

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

# ------------------------------ Convolution Functions -----------------------------
"""
    conv(kG::KGrid, arr1::AbstractVector{ComplexF64}, arr2::AbstractVector{ComplexF64})

Computes the convolution of data over two arrays `arr1` and `arr2`, both given over the fully irreducible BZ.
i.e. ``res[k] = \\sum_{q \\in \\text{BZ}} arr1[k] * arr2[k+q]``.
`arr1` and `arr2` will be expanded to the full BZ internally, the result is given over the fully irreducible BZ.

TODO: at the moment, this function does not use the internal expansion cache of `kG` and is therefore quite slow. 
"""
function conv(kG::KGrid, arr1::AbstractArray{ComplexF64,1}, arr2::AbstractArray{ComplexF64,1})
    Nk(kG) == 1 && return arr1 .* arr2
    tmp = reshape(fft(expandKArr(kG, arr1)) .* fft(expandKArr(kG, arr2)), gridshape(kG)) |> ifft 
    return reduceKArr(kG, ifft_post(kG, tmp)) ./ Nk(kG)
end

"""
    conv!(kG::KGrid, res::AbstractVector{ComplexF64}, arr1::AbstractVector{ComplexF64}, arr2::AbstractVector{ComplexF64})

Inplace version of [`conv`](@ref). The results are written to `res`.
"""
function conv!(kG::KGrid, res::AbstractArray{ComplexF64,1}, arr1::AbstractArray{ComplexF64,1}, arr2::AbstractArray{ComplexF64,1})
    Nk(kG) == 1 && return (res[:] = arr1 .* arr2)
    expandKArr!(kG, arr1)
    tmp = fft(kG.expand_cache)
    expandKArr!(kG, arr2)
    fft!(kG.expand_cache)
    kG.expand_cache[:] = kG.expand_cache .* tmp
    AbstractFFTs.ldiv!(kG.expand_cache, kG.fftw_plan, kG.expand_cache)
    reduceKArr!(kG, res, ifft_post(kG, kG.expand_cache)) 
    res[:] = res ./ Nk(kG)
end


function conv_fft1(kG::KGrid, arr1::AbstractArray{ComplexF64,1}, arr2::AbstractArray{ComplexF64})
    Nk(kG) == 1 && return arr1 .* arr2
    newArr = Array{eltype(arr1),1}(undef, length(kG.kInd))
    conv_fft1!(kG, newArr, arr1, reshape(arr2, gridshape(kG)))
    return newArr
end

function conv_fft1!(kG::KGrid, res::AbstractArray{ComplexF64,1}, arr1::AbstractArray{ComplexF64,1}, arr2::AbstractArray{ComplexF64})
    Nk(kG) == 1 && return (res[:] = arr1 .* arr2)
    expandKArr!(kG, arr1)
    AbstractFFTs.mul!(kG.expand_cache, kG.fftw_plan, kG.expand_cache)
    @simd for i in 1:length(kG.expand_cache)
        @inbounds kG.expand_cache[i] *= arr2[i] 
    end
    AbstractFFTs.ldiv!(kG.expand_cache, kG.fftw_plan, kG.expand_cache)
    reduceKArr!(kG, res, ifft_post(kG, kG.expand_cache)) 
    @simd for i in 1:length(res)
        @inbounds res[i] /= kG.Nk
    end
end

function conv_fft(kG::KGrid, arr1::AbstractArray{ComplexF64}, arr2::AbstractArray{ComplexF64})
    Nk(kG) == 1 && return arr1 .* arr2
    reduceKArr(kG, ifft_post(kG, ifft(arr1 .* arr2))) ./ Nk(kG)
end

function conv_fft!(kG::KGrid, res::AbstractArray{ComplexF64,1}, arr1::AbstractArray{ComplexF64}, arr2::AbstractArray{ComplexF64})
    Nk(kG) == 1 && return (res[:] = arr1 .* arr2)

    @simd for i in eachindex(kG.expand_cache)
        @inbounds kG.expand_cache[i] = arr1[i] .* arr2[i]
    end
    kG.expand_cache[:] = arr1 .* arr2
    AbstractFFTs.ldiv!(kG.expand_cache, kG.fftw_plan, kG.expand_cache)
    reduceKArr!(kG, res, ifft_post(kG, kG.expand_cache)) 
    @simd for i in 1:length(res)
        @inbounds res[i] /= kG.Nk
    end
end
