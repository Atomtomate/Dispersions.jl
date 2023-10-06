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


## Common Functions
The following functions are available for all grids, due to the common struct fields defined above.

TODO

## Interface Functions
