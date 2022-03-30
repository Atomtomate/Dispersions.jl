function gen_kGrid(kg::String, Nk::Int; full=false, sampling=[(2*π/Nk) * j - π for j in 1:Nk])
    sp = findfirst("-", kg)[1]
    data = [kg[1:(sp-1)], kg[(sp+1):end]]
    grid = if lowercase(data[1]) == "3dsc"
        FullKGrid_cP(3, Nk, parse(Float64, data[2]), sampling)
    elseif lowercase(data[1]) == "2dsc"
        FullKGrid_cP(2, Nk, parse(Float64, data[2]), sampling)
    elseif lowercase(data[1]) == "fcc"
        FullKGrid_cF(Nk, parse(Float64, data[2]), [(2*π/Nk) * (j-1) for j in 1:Nk])
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

# ------------------------------ Helper Functions -----------------------------
export reduceKArr, reduceKArr!, expandKArr, expandKArr!, conv, conv!, conv_fft, conv_fft!, conv_fft1, conv_fft1!
"""
    reduceKArr(kGrid::ReducedKGrid, arr)

Takes a kGrid `kGrid` and arbitrary data `arr` over a full BZ and returns an array with data over fully irreducible
BZ. This is mainly used after convolutions, since they require data over the full BZ.
"""
function reduceKArr(kG::ReducedKGrid{gT,D}, arr::AbstractArray{T,D}) where {gT,T,D}
    res = similar(arr, length(kG.kInd))
    reduceKArr!(kG, res, arr)
    return res
end


"""
    reduceKArr!(kGrid::ReducedKGrid, res, arr)

Inplace version of [`reduceKArr`](@ref). The results are written to `res`.
"""
function reduceKArr!(kG::ReducedKGrid{gT,D}, res::AbstractArray{T,1}, arr::AbstractArray{T,D}) where {gT,T,D}
    for (i,ki) in enumerate(kG.kInd)
        @inbounds res[i] = arr[ki...]
    end
end


#TODO: assumed fields: kInd, expand_erms, expand_cache

"""
    expandKArr(kGrid::ReducedKGrid, arr)

Takes a kGrid `kGrid` and arbitrary data `arr` over a reduced BZ and returns an array with data over
full BZ. This is mainly used before convolutions, since they require data over the full BZ.
"""
function expandKArr(kG::ReducedKGrid{gT,D}, arr::AbstractArray{T, 1})::AbstractArray{T, D} where {gT,T,D}
    length(arr) != length(kG.kInd) && throw(ArgumentError("length of k grid ($(length(kG.kInd))) and argument ($(length(arr))) not matching"))
    res = similar(arr, gridshape(kG)...)
    expandKArr!(kG, res, arr)
    return res
end

"""
    expandKArr!(kG, arr)

Inplace version of [`expandKArr`](@ref). The results are written to `kG.expand_cache`.
"""
function expandKArr!(kG::ReducedKGrid, arr::Array{Complex{Float64}, 1})
    for (ri,perms) in enumerate(kG.expand_perms)
        @simd for p in perms
            @inbounds kG.expand_cache[p] = arr[ri]
        end
    end
end

"""
    expandKArr!(kG, [res,] arr)

Inplace version of [`expandKArr`](@ref). The results are written to `res`.
"""
function expandKArr!(kG::ReducedKGrid{gT,D}, res::AbstractArray{T,D}, arr::Array{T, 1}) where {gT,T,D}
    for (ri,perms) in enumerate(kG.expand_perms)
        @simd for p in perms
            @inbounds res[p] = arr[ri]
        end
    end
end


# ------------------------------ Convolution Functions -----------------------------
"""
    conv(kG::ReducedKGrid, arr1::AbstractVector{ComplexF64}, arr2::AbstractVector{ComplexF64})

Computes the convolution of data over two arrays `arr1` and `arr2`, both given over the fully irreducible BZ.
i.e. ``res[k] = \\sum_{q \\in \\text{BZ}} arr1[k] * arr2[k+q]``.
`arr1` and `arr2` will be expanded to the full BZ internally, the result is given over the fully irreducible BZ.

TODO: at the moment, this function does not use the internal expansion cache of `kG` and is therefore quite slow. 
"""
function conv(kG::ReducedKGrid, arr1::AbstractArray{ComplexF64,1}, arr2::AbstractArray{ComplexF64,1})
    Nk(kG) == 1 && return arr1 .* arr2
    tmp = reshape(fft(expandKArr(kG, arr1)) .* fft(expandKArr(kG, arr2)), gridshape(kG)) |> ifft 
    return reduceKArr(kG, ifft_post(kG, tmp)) ./ Nk(kG)
end

"""
    conv!(kG::ReducedKGrid, res::AbstractVector{ComplexF64}, arr1::AbstractVector{ComplexF64}, arr2::AbstractVector{ComplexF64})

Inplace version of [`conv`](@ref). The results are written to `res`.
"""
function conv!(kG::ReducedKGrid, res::AbstractArray{ComplexF64,1}, arr1::AbstractArray{ComplexF64,1}, arr2::AbstractArray{ComplexF64,1})
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


function conv_fft1(kG::ReducedKGrid, arr1::AbstractArray{ComplexF64,1}, arr2::AbstractArray{ComplexF64})
    Nk(kG) == 1 && return arr1 .* arr2
    newArr = Array{eltype(arr1),1}(undef, length(kG.kInd))
    conv_fft1!(kG, newArr, arr1, reshape(arr2, gridshape(kG)))
    return newArr
end

function conv_fft1!(kG::ReducedKGrid, res::AbstractArray{ComplexF64,1}, arr1::AbstractArray{ComplexF64,1}, arr2::AbstractArray{ComplexF64})
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

function conv_fft(kG::ReducedKGrid, arr1::AbstractArray{ComplexF64}, arr2::AbstractArray{ComplexF64})
    Nk(kG) == 1 && return arr1 .* arr2
    reduceKArr(kG, ifft_post(kG, ifft(arr1 .* arr2))) ./ Nk(kG)
end

function conv_fft!(kG::ReducedKGrid, res::AbstractArray{ComplexF64,1}, arr1::AbstractArray{ComplexF64}, arr2::AbstractArray{ComplexF64})
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
