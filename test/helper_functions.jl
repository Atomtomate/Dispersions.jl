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


