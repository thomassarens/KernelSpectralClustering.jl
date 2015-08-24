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
First open the Julia enviroment with a certain number of available processes with _julia -p 24_ (or do _addprocs(23)_ within julia), then specify the used package:
```julia
using KernelSpectralClustering
```
### High Level
Use the __runKSCnetwork__ function to process a network file:
```julia
qTest, baf, k, firstNode, lastNode, trainSubset, numFiles, fileWeighted = runKSCnetwork(networkfile::String, delimiter::Char, mink::Int, maxk::Int)
```
Additional keyword arguments can be supplied:
+ sizeMB=4000::Int -> heuristic to control the size of the binary files and the amount of memory used
+ undirected=true::Bool -> set to false for directed files (extra conversion step using __createMirrored__ function)
+ eval=false::Bool -> set to true to evaluate the quality of the solution using __Coverage__ and __Modularity__ metrics

The function returns the memory-mapped _qTest_ matrix of all nodes in the network and their corresponding community as well as several other variables. The generated binary files are stored in a folder named after the network file in the same directory. A text file _resultKSC.txt_ containing the timings and values of the different steps is also stored in this folder (_metricsKSC.txt_ for optional evaluation).

Consecutive runs on the same dataset can avoid the dataset transformation step by using:
```julia
rerunKSCnetwork(networkfile::String, delimiter::Char, mink::Int, maxk::Int)
```
### Low Level
Use the __kscbignet__ function for networks whose size exceeds available main memory:
```julia
kscbignet(network_file::String, file_dir::String, file_name::String, delimiter::Char, maxMB::Int)
```
Additional keyword arguments can be supplied:
+ fraction=0.15::Float64 -> fraction of network nodes to ne used in subset
+ method="furs"::String -> subset selection method, __randrs__ is a faster, less optimal alternative to __furs__
+ decomp="eigs"::String -> kernel matrix eigen-decomposition method, __eigs__ for approximate, __eig__ for full and __svd__ is planned in future implementation
+ minK=2::Int -> minimum amount of clusters
+ maxK=100::Int -> maximum amount of clusters
+ convertF=true::Bool -> set to false to skip the dataset transdormation step

The __kscnet__ function for in-memory processing of networks is not yet implemented:
```julia
kscnet()
```
