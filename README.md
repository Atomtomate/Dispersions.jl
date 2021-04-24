Types and standard operation for k grids in solid state theory.
This project is inteded to be used as a backend for k grid related operations in order to abstract over the actual k grid logic.

|     Build Status    |      Coverage      |  Documentation |      Social    |
| ------------------- |:------------------:| :-------------:| :-------------:|
| [![Build Status](https://github.com/Atomtomate/Dispersions.jl/workflows/CI/badge.svg)](https://github.com/Atomtomate/Dispersions.jl/actions) | [![Coverage](https://codecov.io/gh/Atomtomate/Dispersions.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/Atomtomate/Dispersions.jl) | [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://Atomtomate.github.io/Dispersions.jl/stable) |[![Gitter](https://badges.gitter.im/JuliansBastelecke/EquivalenceClasses.svg)](https://gitter.im/JuliansBastelecke/EquivalenceClasses?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge) |

# Overview



# Type Hierarchy

In order to provide a common interface, the project defines a rather strict type hierarchy.
The most important supertypes are `FullKGrid{T <: GridType}` and `ReducedKGrid{T <: GridType}` which all `T <: GridType` must implement.
Usually, a certain grid can be defined in different dimensions. These are viewed as separate `GridType` instances and in general do not share functions.

# Interface
