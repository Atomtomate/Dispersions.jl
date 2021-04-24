# Interface

All grid types need to have a certain structure and implement a set of functions which will be discussed here.

## Overview

The interface is divided into common functions which do not need to be implemented for subtypes, since they rely on struct fields all subtypes need to have, and functions which need to be overloaded.

## Common Types

All actual k-grids are supposed to by abstract subtypes of `GridType`. 
For example `abstract type cP_2D <: GridType end` (TODO: ref SC doc). 
Each k-grid of type `T <: GridType` should than implement a `NewKGrid <: FullKGrid{T}` and `NewReducedKGrid <: ReducedKGrid{T} <: ` structs.
Both full and reduced k-grid types are subtypes of `KGrid{T}`.
This hierarchy ensures a set of comon functionalities, that external modules can rely on. 

### Required Fieds
Internally, the structs also need to have at least the following fields:

** Full K-Grid   : **
    - `Nk::Int` that sepcifies the number of k-points, which will be accessed by [`gridPoints`(@ref)
    - `kGrid::T` where `isGridPoints(kG::T) = IsGridPoints{T}` needs to be defined for trait based dispatch to work properly. Generally it should be sufficient to use the predefined `const GridPoints2D = Array{Tuple{Float64,Float64}}` and respective 3D version.
    - `ϵkGrid::Array{Float64, 1}` Array of dispersion relation ``\\epsilon\\_k``.

** Reduced K-Grid: ** 
    - Both fields from the full k-grid need to be present.
    - `kMult::Array{Float64}` specifying the weight of each k point in th reduced grid
    - `kInd::T`  where `isGridInd(kG::T) = IsGridInd{T}` needs to be defined for trait based dispatch to work properly. Generally it should be sufficient to use the predefined `const GridInd2D = Array{Tuple{Int,Int}}` and respective 3D version.


## Common Functions
The following functions are available for all grids, due to the common struct fields defined above.

```@docs
Nk(kG::T) where T <: KGrid
gridPoints(kG::T) where T <: KGrid
```

## Interface Functions
