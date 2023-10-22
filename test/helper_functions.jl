# ---------------------------- Auxilliary Functions -------------------------------
# TODO: avoid reshape
cut_mirror(arr::Base.Iterators.ProductIterator) = cut_mirror(collect(arr))
cut_mirror(arr::Array{T,2}) where {T} =
    arr[ceil(Int, size(arr, 1) / 2):end, ceil(Int, size(arr, 2) / 2):end]
cut_mirror(arr::Array{T,3}) where {T} = arr[
    ceil(Int, size(arr, 1) / 2):end,
    ceil(Int, size(arr, 2) / 2):end,
    ceil(Int, size(arr, 3) / 2):end,
]
#reverse cut. This is a helper function to avoid reversing the array after fft-convolution trick. assumes reserved input and returns correct array including cut
ifft_cut_mirror(arr::Base.Iterators.ProductIterator) = ifft_cut_mirror(collect(arr))
ifft_cut_mirror(arr::Array{T,2}) where {T} =
    arr[end-1:-1:ceil(Int, size(arr, 1) / 2 - 1), end-1:-1:ceil(Int, size(arr, 2) / 2 - 1)]
ifft_cut_mirror(arr::Array{T,3}) where {T} = arr[
    end-1:-1:ceil(Int, size(arr, 1) / 2 - 1),
    end-1:-1:ceil(Int, size(arr, 2) / 2 - 1),
    end-1:-1:ceil(Int, size(arr, 3) / 2 - 1),
]


function test_cut(arr)
    arr2 = copy(arr)
    arr3 = ifft_cut_mirror(arr2)
    ll = floor(Int, size(arr, 1) / 2 + 1)
    index =
        (ndims(arr)) == 2 ? [[x, y] for x = 1:ll for y = 1:x] :
        [[x, y, z] for x = 1:ll for y = 1:x for z = 1:y]
    arr4 = Array{eltype(arr)}(undef, length(index))
    for (i, ti) in enumerate(index)
        arr4[i] = arr3[ti...]
    end
    return arr4
end


function reduce_old(kGrid)
    kGrid_arr = collect(kGrid)
    Nk = size(kGrid_arr, 1)
    if ndims(kGrid_arr) == 2
        index = [[x, y] for x = 1:Nk for y = 1:x]
    elseif ndims(kGrid_arr) == 3
        index = [[x, y, z] for x = 1:Nk for y = 1:x for z = 1:y]
    else
        throw(BoundsError("Number of dimensions for grid must be 2 or 3"))
    end
    grid_red = Array{eltype(kGrid_arr)}(undef, length(index))
    for (i, ti) in enumerate(index)
        grid_red[i] = kGrid_arr[ti...]
    end
    return grid_red
end

# typical use case.
function naive_bubble(fk, pp::Bool = false)
    ωn = 0
    νn = 0
    Σ_loc = 2.734962277113537 - 0.41638191263582125im
    β = 6.0
    μ = 3.512282168483125
    s = pp ? -1 : 1
    disp(ki) = Dispersions.gen_ϵkGrid(T, [ki], fk.t)[1]
    res = zeros(Complex{Float64}, length(fk.kGrid))
    k_en = pp ? reverse(enumerate(fk.kGrid)) : enumerate(fk.kGrid)

    for (kii, ki) in k_en
        for (qii, qi) in enumerate(fk.kGrid)
            Σ_int_ωn = Σ_loc
            Σ_int_ωn_νn = Σ_loc
            w1 = 1im * (π * (2 * ωn + 1) / β) + μ - Σ_int_ωn - disp(ki)
            w2 = 1im * (π * (2 * (ωn + s*νn) + 1) / β) + μ - Σ_int_ωn_νn - disp(qi .+ ki)
            res[qii] = res[qii] - 1.0 / (w1 * w2)
        end
    end
    return res
end

function check_q(vals, Nk; pp=false, transform_back=true)
    s = pp ? +1 : -1
    q_list = Vector{eltype(vals[1])}(undef, length(vals))
    for (i,el) in enumerate(vals)
        q_list[i] =  transform_back ? (mod.(el[2] .+ s .* el[1] .- 2π/Nk .+ π, 2π) .- π .+ 2π/Nk) : (el[2] .+ s .* el[1])
    end
    return q_list
end

function naive_conv(arr1::AbstractArray, arr2::AbstractArray, k0ind::Tuple; pp::Bool=false, round_entries=true)
    res = eltype(arr1) <: Number ? zeros(eltype(arr1), size(arr1)) : nothing
    indices = Array{Vector{Tuple}}(undef, size(arr1))
    res_dbg2 = Array{Vector{Tuple}}(undef, size(arr1))
    arr2_int = pp ? reverse(arr2) : arr2
    s = pp ? -1 : +1
    for i in CartesianIndices(arr1)
        indices[i] = []
        res_dbg2[i] = []
        for j in CartesianIndices(arr1)
            ipj_pre = pp ? Tuple(i) .+ s .* Tuple(j) : Tuple(i) .+ s .* Tuple(j) .- k0ind
            ipj = mod1.(
                    ipj_pre,
                    size(arr1),
                )
            push!(indices[i], (Tuple(j),ipj))
            push!(res_dbg2[i], (arr1[j], arr2[ipj...]))
            eltype(arr1) <: Number && (res[i] += arr1[j] * arr2[ipj...])
        end
    end
    res_dbg2 = round_entries ? map(x->map(xi->map(xii->round.(xii,digits=1),xi),x),res_dbg2) : res_dbg2
    pp && eltype(arr1) <: Number && (res = circshift(res, k0ind))
    pp && (indices = circshift(indices, k0ind))
    pp && (res_dbg2 = circshift(res_dbg2, k0ind))
    # !pp && (res = circshift(res, k0_ind .- 1) )
    # !pp && (indices = circshift(indices, k0_ind .- 1) )
    # !pp && (res_dbg2 = circshift(res_dbg2, k0_ind .- 1) )
    return res, indices, res_dbg2
end

function naive_conv_fft_def(arr1::AbstractArray, arr2::AbstractArray)
    res = zeros(eltype(arr1), size(arr1))
    for j in CartesianIndices(arr1)
        for i in CartesianIndices(arr2)
            ii = mod.(Tuple(i) .- (Tuple(j)), size(arr2)) .+ Tuple(ones(Int, length(i)))
            res[i] += arr1[j] * arr2[ii...]
        end
    end
    return res
end

function fcc_dispersion(k::Tuple, t::Float64)
    return -4t*(cos(k[1]/2) * cos(k[2]/2) + cos(k[1]/2) * cos(k[3]/2) + cos(k[2]/2) * cos(k[3]/2))
end

transform_k_to_minus_k(k, kG) = map(ki -> ki != 0 ? mod.(-1 .* ki .+ π .- 1/Nk(kG), 2π) .- π .+ 1/Nk(kG) : ki, k)
