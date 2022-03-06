function gen_kGrid(kg::String, Nk::Int; rot_angles=nothing)
    #TODO: check for reasonable input
    sp = findfirst("-", kg)[1]
    data = [kg[1:(sp-1)], kg[(sp+1):end]]
    grid = if lowercase(data[1]) == "3dsc"
        KGrid(SC, 3, Nk, parse(Float64, data[2]), rot_angles=rot_angles)
    elseif lowercase(data[1]) == "2dsc"
        KGrid(SC, 2, Nk, parse(Float64, data[2]), rot_angles=rot_angles)
    elseif lowercase(data[1]) == "fcc"
        KGrid(FCC, 3, Nk, parse(Float64, data[2]), rot_angles=rot_angles)
    elseif lowercase(data[1]) == "p6m"
        KGrid_p6m(p6m, 2, Nk, parse(Float64, data[2]), rot_angles=rot_angles)
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

Computes the convolution with a plus sign of data over two arrays `arr1` and `arr2`,
i.e. ``res[k] = \\sum_{q \\in \\text{BZ}} arr1[k] * arr2[k+q]``.
"""
function conv(kG::KGrid, arr1::AbstractArray{ComplexF64,1}, arr2::AbstractArray{ComplexF64,1})
    Nk(kG) == 1 && return arr1 .* arr2
    res = similar(arr2)
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
    AbstractFFTs.mul!(kG.fft_cache, kG.fftw_plan, reverse(reshape(view(arr2,:),gs)))
    conv_fft!(kG, res, res, kG.fft_cache)
end

"""
    conv_fft1(kG::KGrid, res::AbstractVector{ComplexF64}, arr1::AbstractVector{ComplexF64}, arr2::AbstractVector{ComplexF64})

Version of [`conv`](@ref) for precomputed `arr2`. Note, that reversing the array is also
expected, i.e. `arr2 = fft(reverse(in_arr2))`.
"""
function conv_fft1(kG::KGrid, arr1::AbstractArray{ComplexF64}, arr2::AbstractArray{ComplexF64})
    Nk(kG) == 1 && return arr1 .* arr2
    res = similar(arr2)
    conv_fft1!(kG, res, arr1, arr2)
    return res
end

"""
    conv_fft1!(kG::KGrid, res::AbstractVector{ComplexF64}, arr1::AbstractVector{ComplexF64}, arr2::AbstractVector{ComplexF64})

Inplace version of [`conv_fft1`](@ref).
"""
function conv_fft1!(kG::KGrid, res::AbstractArray{ComplexF64,1}, arr1::AbstractArray{ComplexF64,1}, arr2::AbstractArray{ComplexF64})
    Nk(kG) == 1 && return (res[:] = arr1 .* arr2)
    AbstractFFTs.mul!(kG.fft_cache, kG.fftw_plan, reshape(view(arr1,:),gridshape(kG)))
    conv_fft!(kG, view(res,:), view(kG.fft_cache,:), view(arr2,:))
end

"""
    conv_fft(kG::KGrid, res::AbstractVector{ComplexF64}, arr1::AbstractVector{ComplexF64}, arr2::AbstractVector{ComplexF64})

Version of [`conv`](@ref) for precomputed `arr1` and `arr2`. Note, that reversing `arr2` is also
expected, see also [`conv_fft1`](@ref).
"""
function conv_fft(kG::KGrid, arr1::AbstractArray{ComplexF64}, arr2::AbstractArray{ComplexF64})
    Nk(kG) == 1 && return arr1 .* arr2
    res = similar(arr2)
    conv_fft!(kG, view(res,:), view(arr1,:), view(arr2,:))
    return res
end


"""
    conv_fft!(kG::KGrid, res::AbstractVector{ComplexF64}, arr1::AbstractVector{ComplexF64}, arr2::AbstractVector{ComplexF64})

Inplace version of [`conv_fft`](@ref).
"""
function conv_fft!(kG::KGrid, res::AbstractArray{ComplexF64,1}, arr1::AbstractArray{ComplexF64}, arr2::AbstractArray{ComplexF64})
    Nk(kG) == 1 && return (res[:] = arr1 .* arr2)
    res_v = reshape(view(res,:),gridshape(kG))
    for i in eachindex(kG.fft_cache)
        @inbounds kG.fft_cache[i] = arr1[i] .* arr2[i]
    end
    AbstractFFTs.ldiv!(res_v, kG.fftw_plan, kG.fft_cache)
    ifft_post!(kG, res_v)
    res[:] = res_v[:] ./ Nk(kG)
end
