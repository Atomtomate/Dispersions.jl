# ---------------------------- Auxilliary Functions -------------------------------
# TODO: avoid reshape
cut_mirror(arr::Base.Iterators.ProductIterator) = cut_mirror(collect(arr))
cut_mirror(arr::Array{T, 2}) where T = arr[ceil(Int,size(arr,1)/2):end, ceil(Int,size(arr,2)/2):end]
cut_mirror(arr::Array{T, 3}) where T = arr[ceil(Int,size(arr,1)/2):end, ceil(Int,size(arr,2)/2):end, ceil(Int,size(arr,3)/2):end]
#reverse cut. This is a helper function to avoid reversing the array after fft-convolution trick. assumes reserved input and returns correct array including cut
ifft_cut_mirror(arr::Base.Iterators.ProductIterator) = ifft_cut_mirror(collect(arr))
ifft_cut_mirror(arr::Array{T, 2}) where T = arr[end-1:-1:ceil(Int,size(arr,1)/2-1), 
                                                          end-1:-1:ceil(Int,size(arr,2)/2-1)]
ifft_cut_mirror(arr::Array{T, 3}) where T = arr[end-1:-1:ceil(Int,size(arr,1)/2-1), 
                                                          end-1:-1:ceil(Int,size(arr,2)/2-1), 
                                                          end-1:-1:ceil(Int,size(arr,3)/2-1)]


function test_cut(arr)
    arr2 = copy(arr)
    arr3 = ifft_cut_mirror(arr2)
    ll = floor(Int,size(arr,1)/2 + 1)
    index = (ndims(arr)) == 2 ? [[x,y] for x=1:ll for y=1:x] : [[x,y,z] for x=1:ll for y=1:x for z = 1:y]
    arr4 = Array{eltype(arr)}(undef, length(index))
    for (i,ti) in enumerate(index)
        arr4[i] = arr3[ti...]
    end
    return arr4
end


function reduce_old(kGrid) 
    kGrid_arr = collect(kGrid)
    Nk = size(kGrid_arr,1)
    if ndims(kGrid_arr) == 2
        index = [[x,y] for x=1:Nk for y=1:x]
    elseif ndims(kGrid_arr) == 3
        index = [[x,y,z] for x=1:Nk for y=1:x for z = 1:y]
    else
        throw(BoundsError("Number of dimensions for grid must be 2 or 3"))
    end
    grid_red = Array{eltype(kGrid_arr)}(undef, length(index))
    for (i,ti) in enumerate(index)
        grid_red[i] = kGrid_arr[ti...]
    end
    return grid_red
end

# typical use case.
function naive_bubble(fk::KGrid{T}) where T
    ωn = 0
    νn = 0
    Σ_loc = 2.734962277113537 - 0.41638191263582125im
    β = 6.0
    μ = 3.512282168483125
    disp(ki) = Dispersions.gen_ϵkGrid(T, [ki], fk.t)[1]
    res = zeros(Complex{Float64}, length(fk.kGrid))

    for (kii,ki) in enumerate(fk.kGrid)
        for (qii,qi) in enumerate(fk.kGrid)
            Σ_int_ωn = Σ_loc
            Σ_int_ωn_νn = Σ_loc
            w1 = 1im * (π*(2*ωn+1)/β) + μ - Σ_int_ωn - disp(ki)
            w2 = 1im * (π*(2*(ωn+νn)+1)/β) + μ - Σ_int_ωn_νn - disp(ki .+ qi)
            res[qii] = res[qii] - 1.0 / (w1 * w2)
        end
    end
    return res
end

# function naive_conv(kG1::ReducedKGrid, kG2::ReducedKGrid, arr1::AbstractArray, arr2::AbstractArray)
#     Nk(kG1) == 1 && return arr1 .* arr2
#     a1 = expandKArr(kG1, arr1)
#     a2 = expandKArr(kG2, arr2)
#     res = zeros(eltype(arr1), size(a1))
#     for i in CartesianIndices(a1)
#         for j in CartesianIndices(a2)
#             ii = mod1.(Tuple(j) .- Tuple(i) .- Tuple(ones(Int,length(i))), size(a1))
#             println("res[$(i.I)] += a1[$(j.I)] * a2[$(ii)] = $(round.(a1[i]*a2[ii...],digits=3))")
#             res[i] += a1[i]*a2[ii...]
#         end
#     end
#     return res ./ Nk(kG1)
# end

# function naive_conv(kG::ReducedKGrid, arr1::AbstractArray, arr2::AbstractArray)
#     Nk(kG) == 1 && return arr1 .* arr2
#     a1 = expandKArr(kG, arr1)
#     a2 = expandKArr(kG, arr2)
#     res = zeros(eltype(arr1), size(a1))
#     for j in CartesianIndices(a2)
#         for i in CartesianIndices(a1)
#             ii = mod1.(Tuple(j) .+ Tuple(i) .- Tuple(ones(Int,length(i))), size(a1))
#             res[j] += a1[j]*a2[ii...]
#         end
#     end
#     return res ./ Nk(kG)
# end
