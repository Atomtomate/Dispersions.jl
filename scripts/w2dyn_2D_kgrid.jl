using Pkg
Pkg.activate("/scratch/projects/hhp00048/codes/Dispersions.jl")

using Dispersions
kgrid_str = ARGS[1]
N    = parse(Int, ARGS[2])
outf = ARGS[3]

#kG = gen_kGrid("2Dsc-0.25-0.075-0.05", N)
kG = gen_kGrid(kgrid_str, N)
d = expandKArr(kG,dispersion(kG))[:]
k = expandKArr(kG,kG.kGrid)[:]

open(outf, "w") do f
    write(f, "$(length(d))\t1\t1\n")
    for i in 1:length(d)
        write(f, "$(round(k[i][1]/(2*π),digits=6))\t$(round(k[i][2]/(2*π),digits=6))\t0.0\n")
        write(f, "$(round(d[i],digits=6))\t0.0\n")
    end
end
