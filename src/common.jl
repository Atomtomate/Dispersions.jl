function gen_kGrid(kg::String, Nk::Int)
    data = split(kg, "-")
    grid = if lowercase(data[1]) == "3dsc"
        FullKGrid_cP_3D(Nk, parse(Float64, data[2]))
    elseif lowercase(data[1]) == "2dsc"
        FullKGrid_cP_2D(Nk, parse(Float64, data[2]))
    else
        throw(ArgumentError("Unkown grid type"))
    end
    return grid
end
