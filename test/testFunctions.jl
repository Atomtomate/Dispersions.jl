function naive_conv(kG::ReducedKGrid, arr1::AbstractArray{T,1}, arr2::AbstractArray{T,1}) where T <: Number
    Nk(kG) == 1 && return arr1 .* arr2
    a1 = expandKArr(kG, arr1)
    a2 = expandKArr(kG, arr2)
    res = zeros(size(a1)...)
    mi = size(a2)
    for j in CartesianIndices(a2)
        for i in CartesianIndices(a1)
            ii = mod1.(Tuple(j) .+ Tuple(i) .- Tuple(ones(Int,length(i))), size(a1))
            res[j] += a1[j]*a2[ii...]
            println("j=$j, i=$i", a1[j], " * ",a2[ii...], " = ", res[j])
        end
    end
    return res
end
function naive_conv(kG::ReducedKGrid, arr1::AbstractArray{T,2}, arr2::AbstractArray{T,2}) where T <: Number
    Nk(kG) == 1 && return arr1 .* arr2
    a1 = arr1
    a2 = arr2
    res = zeros(size(a1)...)
    mi = size(a2)
    for j in CartesianIndices(a2)
        for i in CartesianIndices(a1)
            ii = mod1.(Tuple(j) .+ Tuple(i) .- Tuple(ones(Int,length(i))), size(a1))
            res[j] += a1[j]*a2[ii...]
            println("j=$j, i=$i", a1[j], " * ",a2[ii...], " = ", res[j])
        end
    end
    return res
end
