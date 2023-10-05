Types and standard operation for k grids and tight binding dispersions in solid state theory.
This project is inteded to be used as a backend for k grid related operations in order to abstract over the actual k grid logic.
For an example usage see [LadderDGA.jl](https://github.com/Atomtomate/LadderDGA.jl)

|     Build Status    |      Coverage      |  Documentation |      Social    |
| ------------------- |:------------------:| :-------------:| :-------------:|
| [![Build Status](https://github.com/Atomtomate/Dispersions.jl/workflows/CI/badge.svg)](https://github.com/Atomtomate/Dispersions.jl/actions) | [![codecov](https://codecov.io/gh/Atomtomate/Dispersions.jl/branch/master/graph/badge.svg?token=FbJKjHb7DW)](https://codecov.io/gh/Atomtomate/Dispersions.jl) | [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://Atomtomate.github.io/Dispersions.jl/stable) |[![Gitter](https://badges.gitter.im/JuliansBastelecke/EquivalenceClasses.svg)](https://gitter.im/JuliansBastelecke/EquivalenceClasses?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge) |

# Overview


# Interface

- `kG = genKGrid(s, nk)` will generate a k grid with `nk` points in each direction from a configuration string. Examples are `s = "2Dsc-1.0"` (2D simple cubic with hopping parameter - `t = 1.0`), `s = "FCC-1.0-1.1-1.3"` for a FCC lattice.
- `gridshape(kG)` returns the shape of a k grid. Example: `(4,4)` for a lattice generated with `genKGrid("2Dsc-1.0", 4)`
- `Nk(kG)` returns the number of k points for the given grid. Example: `16` for a lattice generated with `genKGrid("2Dsc-1.0", 4)`
- `gridPoints(kG)` returns the k vektors for the given grid.
- `dispersion(kG)` returns the dispersion relation, i.e. an array of energy value for each point of `gridPoints(kG)`.
- `grid_type(kG)`  returns the KGridType of the given grid without the dimension the grid lives on.

- `conv(kG, arr1, arr2)` returns the convolution of the arrays `arr1` and `arr2` (any function evaluated on `kG`) on `kG`. This is usefull for the evaluation of expressions like G(q) = sum_k G(q)*G(k+q)
- `conv_fft(kG, arr1, arr2)` and `conv_fft1(kG, arr1, arr2)` are helpers, used when `arr1` or both `arr1` and `arr2` are already fourier transformed (used for repeated evaluation with the same kernel).
