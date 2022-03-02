function gen_kGrid(kg::String, Nk::Int)
    #TODO: check for reasonable input
    sp = findfirst("-", kg)[1]
    data = [kg[1:(sp-1)], kg[(sp+1):end]]
    grid = if lowercase(data[1]) == "3dsc"
        KGrid(SC, 3, Nk, parse(Float64, data[2]))
    elseif lowercase(data[1]) == "2dsc"
        KGrid(SC, 2, Nk, parse(Float64, data[2]))
    elseif lowercase(data[1]) == "fcc"
        KGrid(FCC, 3, Nk, parse(Float64, data[2]))
    elseif lowercase(data[1]) == "p6m"
        KGrid_p6m(p6m, 2, Nk, parse(Float64, data[2]))
    else
        throw(ArgumentError("Unkown grid type"))
    end
    return grid
end



# ================================================================================ #
#                                   Interface                                      #
# ================================================================================ #
export conv, conv!, conv_fft, conv_fft!, conv_fft1, conv_fft1!

# ------------------------------ Helper Functions -----------------------------

"""
    gridshape(kG::T) where T <: KGrid

shape of kGrid (e.g. `(kG.Ns, kG.Ns)` for 2D sc) 
"""
gridshape(kG::KGrid{T,D}) where {T,D} = ntuple(_ -> kG.Ns, D)

"""
    Nk(kG::T) where T <: KGrid

Total number of k points.
"""
Nk(kG) = kG.Nk

"""
    gridPoints(kG::T)::Int where T <: KGrid

k vectors of grid.
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
    res = similar(arr1)
    conv!(kG, res, arr1, arr2)
    return res
end

"""
    conv!(kG::KGrid, res::AbstractVector{ComplexF64}, arr1::AbstractVector{ComplexF64}, arr2::AbstractVector{ComplexF64})

Inplace version of [`conv`](@ref). The results are written to `res`.
"""
function conv!(kG::KGrid, res::AbstractArray{ComplexF64,1}, arr1::AbstractArray{ComplexF64,1}, arr2::AbstractArray{ComplexF64,1})
    Nk(kG) == 1 && return (res[:] = arr1 .* arr2)
    gs = gridshape(kG)
    res_v = reshape(view(res,:),gs)
    AbstractFFTs.mul!(res_v, kG.fftw_plan, reshape(view(arr1,:),gs))
    AbstractFFTs.mul!(kG.fft_cache, kG.fftw_plan, reshape(view(arr2,:),gs))
    @inbounds res_v = res_v .* kG.fft_cache
    AbstractFFTs.ldiv!(kG.fft_cache, kG.fftw_plan, res_v)
    ifft_post!(kG, res_v, kG.fft_cache)
    res[:] = res_v[:] ./ Nk(kG)
end


function conv_fft1(kG::KGrid, arr1::AbstractArray{ComplexF64}, arr2::AbstractArray{ComplexF64})
    Nk(kG) == 1 && return arr1 .* arr2
    res = Array{eltype(arr1),1}(undef, length(arr1))
    conv_fft1!(kG, res, arr1, arr2)
    return res
end

function conv_fft1!(kG::KGrid, res::AbstractArray{ComplexF64,1}, arr1::AbstractArray{ComplexF64,1}, arr2::AbstractArray{ComplexF64})
    Nk(kG) == 1 && return (res[:] = arr1 .* arr2)
    AbstractFFTs.mul!(kG.fft_cache, kG.fftw_plan, reshape(view(arr1,:),gridshape(kG)))
    conv_fft!(kG, view(res,:), view(kG.fft_cache,:), view(arr2,:))
end

function conv_fft(kG::KGrid, arr1::AbstractArray{ComplexF64}, arr2::AbstractArray{ComplexF64})
    Nk(kG) == 1 && return arr1 .* arr2
    res = similar(arr1)
    conv_fft!(kG, view(res,:), view(arr1,:), view(arr2,:))
    return res
end

function conv_fft!(kG::KGrid, res::AbstractArray{ComplexF64,1}, arr1::AbstractArray{ComplexF64}, arr2::AbstractArray{ComplexF64})
    Nk(kG) == 1 && return (res[:] = arr1 .* arr2)
    res_v = reshape(view(res,:),gridshape(kG))
    @simd for i in eachindex(kG.fft_cache)
        @inbounds res[i] = arr1[i] .* arr2[i]
    end
    AbstractFFTs.ldiv!(kG.fft_cache, kG.fftw_plan, res_v)
    ifft_post!(kG, res_v, kG.fft_cache)
    res[:] = res_v[:] ./ Nk(kG)
    # @simd for i in eachindex(res_v)
    #     @inbounds res_v[i] /= kG.Nk
    # end
end
