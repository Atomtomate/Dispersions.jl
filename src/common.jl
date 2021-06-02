function gen_kGrid(kg::String, Nk::Int)
    data = split(kg, "-")
    grid = if lowercase(data[1]) == "3dsc"
        FullKGrid_cP_3D(Nk, parse(Float64, data[2]))
    elseif lowercase(data[1]) == "2dsc"
        FullKGrid_cP_2D(Nk, parse(Float64, data[2]))
    elseif lowercase(data[1]) == "p6m"
        FullKGrid_p6m(Nk, parse(Float64, data[2]))
    else
        throw(ArgumentError("Unkown grid type"))
    end
    return grid
end

function conv(kG::ReducedKGrid, arr1::AbstractArray{Complex{Float64},1}, arr2::AbstractArray{Complex{Float64},1})
    kG.Nk == 1 && return arr1 .* arr2
    reshape(fft(expandKArr(kG, arr1)) .* fft(expandKArr(kG, arr2)), gridshape(kG)) |> ifft |> x-> reduceKArr_reverse(kG, x) ./ kG.Nk
end


function conv_fft1(kG::ReducedKGrid, arr1::AbstractArray{Complex{Float64},1}, arr2::AbstractArray{Complex{Float64},1})
    kG.Nk == 1 && return arr1 .* arr2
    reshape(fft(expandKArr(kG, arr1))[:] .* arr2, gridshape(kG)) |> ifft |> x-> reduceKArr_reverse(kG, x) ./ kG.Nk
end

function conv_fft(kG::ReducedKGrid, arr1::AbstractArray{Complex{Float64},1}, arr2::AbstractArray{Complex{Float64},1})
    kG.Nk == 1 && return arr1 .* arr2
    reshape(arr1 .* arr2, gridshape(kG)...) |> ifft |> x-> reduceKArr_reverse(kG, x) ./ kG.Nk
end
