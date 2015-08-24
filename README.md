# KernelSpectralClustering

[![Build Status](https://travis-ci.org/thomassarens/KernelSpectralClustering.jl.svg?branch=master)](https://travis-ci.org/thomassarens/KernelSpectralClustering.jl)

## Installation
Use the Julia package manager:
```julia
Pkg.clone("git://github.com/thomassarens/KernelSpectralClustering.jl.git")
```
## Dataset Format
The __KernelSpectralClustering__ module is designed to apply community detection to un-/directed networks formated as textfiles with 2 (unweighted) or 3 (weighted) delimiter sparated values.

## Usage
First open the Julia enviroment with a certain number of available processes with 'julia -p 24' (or do 'addprocs(23)' within julia), then specify the used package:
```julia
using KernelSpectralClustering
```
```julia
runKSCnetwork(networkfile::String, delimiter::Char, mink::Int, maxk::Int)
```
Additional keyword arguments can be used:
*sizeMB=4000::Int -> heuristic to control the size of the binary files and the amount of memory used
*undirected=true::Bool -> set to false for directed files (extra conversion step using __createMirrored__ function)
*eval=false::Bool -> set to true to evaluate the quality of the solution using __Coverage__ and __Modularity__ metrics

Consecutive runs on the same dataset can avoid the dataset transformation step by using:
```julia
rerunKSCnetwork(networkfile::String, delimiter::Char, mink::Int, maxk::Int)
```
