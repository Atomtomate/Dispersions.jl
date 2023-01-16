# ================================================================================ #
#                                     Hexagonal                                    #
# ================================================================================ #

abstract type p6m <: KGridType end

# -------------------------------------------------------------------------------- #
#                            Grid Generator Functions                              #
# -------------------------------------------------------------------------------- #

gen_sampling(::Type{p6m}, D::Int, Ns::Int) = Base.product([[(2*j - Ns - 1)/(2*Ns) for j in 1:Ns] for Di in 1:2]...)
basis_transform(::Type{p6m}, v::Tuple) = (2π * v[1], -2π * v[1]/sqrt(3) + 4π*v[2]/sqrt(3))

function reduce_KGrid(::Type{p6m}, D::Int, Ns::Int, kGrid::AbstractArray)
    ind = collect(Base.product([1:Ns for Di in 1:2]...))
    kMult = ones(length(ind))
    expand_perms = map(x -> [CartesianIndex{2}(x)],ind[:])
    red_map = CartesianIndex.(ind[:])
    red_conv_map = reverse(red_map)
    return red_map, red_conv_map, kMult, expand_perms, kGrid[:]
end

function gen_ϵkGrid(::Type{p6m}, kGrid::GridPoints, t::T, tp::T, tpp::T) where {T<:Real}
    if tp != 0.0 || tpp != 0.0
        throw(ArgumentError("Dispersion of p6m not implemented for non-zero next nearest neightbor hopping!"))
    end
    gen_ϵkGrid(p6m, kGrid, t) 
end

gen_ϵkGrid(::Type{p6m}, kGrid::GridPoints{2}, t::T1) where T1 <: Number = collect(map(kᵢ -> -2*t*(cos.(0.5*(kᵢ[1] + sqrt(3)*kᵢ[2])) + cos(0.5*(kᵢ[1] - sqrt(3)*kᵢ[2])) + cos(kᵢ[1])), kGrid))

# -------------------------------------------------------------------------------- #
#                             Custom Helper Functions                              #
# -------------------------------------------------------------------------------- #

conv_sample_post(kG::KGrid{p6m, 2}, x) = x
