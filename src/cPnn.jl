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

function reduce_KGrid(::Type{cPnn}, D::Int, Ns::Int, kGrid::AbstractArray)
    ind = collect(Base.product([1:Ns for Di in 1:2]...))
    kMult = ones(length(ind))
    expand_perms = map(x -> [CartesianIndex{2}(x)],ind[:])
    red_map = CartesianIndex.(ind[:])
    red_conv_map = reverse(red_map)
    return red_map, red_conv_map, kMult, expand_perms, kGrid[:]
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
conv_sample_post(kG::KGrid{cPnn, 2}, x) = x
