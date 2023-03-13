# ================================================================================ #
#                            2D Cubic for Hofstadter Model                         #
# ================================================================================ #

abstract type Hofstadter{P,Q} <: KGridType end

# -------------------------------------------------------------------------------- #
#                            Grid Generator Functions                              #
# -------------------------------------------------------------------------------- #

gen_sampling(::Type{Hofstadter{P,Q}}, D::Int, Ns::Int) where {P,Q} = Base.product([[(2*j - Ns)/(2*Ns) for j in 0:Ns-1] for Di in 1:2]...)
basis_transform(::Type{Hofstadter{P,Q}}, v::Tuple) where {P,Q} = (2π * v[1], 2π*v[2]/Q)

function reduce_KGrid(::Type{Hofstadter{P,Q}}, D::Int, Ns::Int, kGrid::AbstractArray) where {P,Q}
    ind = collect(Base.product([1:Ns for Di in 1:2]...))
    kMult = ones(length(ind))
    expand_perms = map(x -> [CartesianIndex{2}(x)],ind[:])
    red_map = CartesianIndex.(ind[:])
    red_conv_map = reverse(red_map)
    return red_map, red_conv_map, kMult, expand_perms, kGrid[:]
end

function gen_ϵkGrid_Hofstadter_ij(::Type{Hofstadter{P,Q}}, i::Int, j::Int, kᵢ::NTuple{2,Float64}, t::Float64) where {P,Q}
    if i==j
        return -2*t*(cos(kᵢ[1] + 2*π*P*(i-1)/Q))
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

function gen_ϵkGrid(::Type{Hofstadter{P,Q}}, kGrid::GridPoints, t::Float64, tp::Float64, tpp::Float64) where {P,Q}
    if gcd(P,Q) != 1 || Q < 3
        throw(ArgumentError("gcd(P,Q) = $(gcd(P,Q)) != 1 or Q = $(Q) < 3!"))
    end
    res = Array{ComplexF64,3}(undef, Q, Q, length(kGrid))
    for (l, kᵢ) in enumerate(kGrid)
        for j in 1:Q
            for i in 1:Q
                res[i,j,l] = gen_ϵkGrid_Hofstadter_ij(Hofstadter{P,Q}, i, j, kᵢ, t)
            end
        end
    end
    res
end

# -------------------------------------------------------------------------------- #
#                             Custom Helper Functions                              #
# -------------------------------------------------------------------------------- #

conv_sample_post(kG::KGrid{Hofstadter{P,Q}, 2}, x) where {P,Q} = x
