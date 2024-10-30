# ================================================================================ #
#                                   Simple Cubic                                   #
# ================================================================================ #

abstract type cPnn <: KGridType end

# -------------------------------------------------------------------------------- #
#                            Grid Generator Functions                              #
# -------------------------------------------------------------------------------- #

gen_sampling(::Type{cPnn}, D::Int, Ns::Int) =
    Base.product([[(2 * π / Ns) * j - π for j = 1:Ns] for Di = 1:D]...)

basis_transform(::Type{cPnn}, v::Tuple) = v

transform_to_first_BZ(kG::KGrid{cPnn,D}, k) where D =
    Tuple(map(ki -> ki != 0 ? mod(ki + (π - 1/kG.Nk), 2π) + (-π) + 1/kG.Nk : ki, k))

conv_Indices(::Type{cPnn}, D::Int, Ns::Int) = conv_Indices(cP, D, Ns)

function reduce_KGrid(::Type{cPnn}, D::Int, Ns::Int, kGrid::AbstractArray)
    return reduce_KGrid(cP, D::Int, Ns::Int, kGrid::AbstractArray)
end

# 10.1103/PhysRevLett.87.047003
gen_ϵkGrid(::Type{cPnn}, kGrid::GridPoints, t::T, tp::T, tpp::T) where {T<:Real} =
collect(map(kᵢ -> -2 * t   * sum(cos.(kᵢ)) +
                  4 * tp  * cos(kᵢ[1])*cos(kᵢ[2]) - # TODO: np.einsum like impl. for 3D
                  2 * tpp * sum(cos.(2 .* kᵢ))
            , kGrid))

# -------------------------------------------------------------------------------- #
#                             Custom Helper Functions                              #
# -------------------------------------------------------------------------------- #
conv_sample_post(kG::KGrid{cPnn, 2}, x) =
    ShiftedArrays.circshift(x, floor.(Int, gridshape(kG) ./ 2) .- 1)
#TODO: this somehow works when not doing the reverse on the second input. We should find out why, this makes the convolution a lot faster
conv_post_old(kG::KGrid{cPnn, 2}, x::Array{T,D}) where {D,T<:Number} =
    reduceKArr(kG, ShiftedArrays.circshift(x, floor.(Int, gridshape(kG) ./ 2) .+ 1))
