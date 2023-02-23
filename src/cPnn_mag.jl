# ================================================================================ #
#                            2D Cubic for Hofstadter Model                         #
# ================================================================================ #

abstract type Hofstadter{Q} <: KGridType end

# -------------------------------------------------------------------------------- #
#                            Grid Generator Functions                              #
# -------------------------------------------------------------------------------- #

gen_sampling(::Type{Hofstadter}, D::Int, Ns::Int) = Base.product([[(2*j - Ns)/(2*Ns) for j in 0:Ns-1] for Di in 1:2]...)
basis_transform(::Type{Hofstadter{Q}}, v::Tuple) where Q = (2π * v[1], 2π*v[2]/Q)

function reduce_KGrid(::Type{Hofstadter}, D::Int, Ns::Int, kGrid::AbstractArray)
    ind = collect(Base.product([1:Ns for Di in 1:2]...))
    kMult = ones(length(ind))
    expand_perms = map(x -> [CartesianIndex{2}(x)],ind[:])
    red_map = CartesianIndex.(ind[:])
    red_conv_map = reverse(red_map)
    return red_map, red_conv_map, kMult, expand_perms, kGrid[:]
end

function gen_ϵkGrid_Hofstadter_ij(::Type{Hofstadter{Q}}, i::Int, j::Int, kᵢ::NTuple{2,Float64}, t::Float64) where Q
    if i==j
        return -2*t*(cos(kᵢ[1] + 2*π*(i-1)/Q))
    elseif i == j+1 || i == j-1
        return -t
    elseif i == Q && j == 1
        return t*exp(-1im*Q*kᵢ[2])
    elseif j == Q && i == 1
        return t*exp(1im*Q*kᵢ[2])
    else
        return 0
    end
end

function gen_ϵkGrid(::Type{Hofstadter}, kGrid::GridPoints, t::Float64, p::Int, q::Int)
    if gcd(p,q) != 1 || q < 3
        throw(ArgumentError("gcd(p,q) = $(gcd(p,q)) != 1 or q = $(q) < 3!"))
    end

    gen_ϵkGrid(Hofstadter, kGrid, t) 
end

gen_ϵkGrid(::Type{Hofstadter}, kGrid::GridPoints{2}, t::T1) where T1 <: Number = collect(map(kᵢ -> -2*t*(cos.(0.5*(kᵢ[1] + sqrt(3)*kᵢ[2])) + cos(0.5*(kᵢ[1] - sqrt(3)*kᵢ[2])) + cos(kᵢ[1])), kGrid))

# -------------------------------------------------------------------------------- #
#                             Custom Helper Functions                              #
# -------------------------------------------------------------------------------- #

conv_sample_post(kG::KGrid{Hofstadter, 2}, x) = x
