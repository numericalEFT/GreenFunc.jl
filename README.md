# GreenFunc



[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://numericalEFT.github.io/GreenFunc.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://numericalEFT.github.io/GreenFunc.jl/dev)
[![Build Status](https://github.com/numericalEFT/GreenFunc.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/numericalEFT/GreenFunc.jl/actions/workflows/CI.yml?query=branch%3Amaster)
[![Coverage](https://codecov.io/gh/numericalEFT/GreenFunc.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/numericalEFT/GreenFunc.jl)

GreenFunc.jl provides a general container for multidimensional Green's functions. 
## Features
 - Universal structure for data array that is saved on meshes
 - Structured temporal meshes for imaginary time/Matsubara frequency space.
 - Discrete Lehmann representaion (DLR) meshes, which is a generic and compact representation of Green's functions proposed in the Ref. [1]. Allows fast and accurate Fourier transform between the imaginary-time domain and the Matsubara-frequency domain.
 - Structure of direct product of arbitrary number of meshes.
 - Conversion between Green's function in this package and in TRIQS.
## Installation
This package has been registered. So, simply type `import Pkg; Pkg.add("GreenFunc")` in the Julia REPL to install.

## Basic Usage
We first show how to generate a multidimensional function given by the following equation: $G_{\sigma_1,\sigma_2}(q,\omega_n) = \frac{1}{i \omega_n - \epsilon}$. This is a bare green's function of particles that have two inner spin 1/2 degrees, $sigma_1$ and $sigma_2$, and a dispersion $\epsilon = q^2-1$.
     
```julia
julia> using GreenFunc
julia> using CompositeGrids
julia> mesh1 = SimpleGrid.Uniform{Float64}([0.0, 10.0], 50);
julia> mesh2 = MeshGrids.ImFreq(100.0, FERMION);
julia> innermesh = (1:2,1:2);    
julia> g =  MeshArray(innermesh..., mesh1, mesh2, dtype=ComplexF64)
Meshed array with dims = (2, 2, 50, 60) and total length = 12000
- Mesh: Tuple{UnitRange{Int64}, UnitRange{Int64}, CompositeGrids.SimpleG.Uniform{Float64}, ImFreq{Float64, CompositeGrids.SimpleG.Arbitrary{Int64}}}
julia> for ind in eachindex(g)
           q = g.mesh[3][ind[3]]
           ω_n = g.mesh[4][ind[4]]
           g[ind] = 1/(im*ω_n - (q^2-1))
       end
```
- Momentum grids are handled by CompositeGrids package. Here SimpleGrid.Uniform{T}([mink,maxk],N) generates a linearly spaced momentum grid with N points in the interval (mink,maxk). For other specialized grids see https://github.com/numericalEFT/CompositeGrids.jl.
- GreenFunc package provides three specialized grids: imaginary time, matsubara frequnecy and DLR grids. They can be generated with ImTime, ImFreq and DLRFreq(beta,statistics) methods correspondingly. Here beta is the inverse temperature, and statistics defines whether the particle is fermion or boson. 
- By default the data type of MeshArray is Float64, and data is set to nothing if not provided. One can reset them with 'dtype' and 'data' keywards arguments. 
- Meshes for all variables, including the discrete ones (spin) and continuous ones (momentum and frequnecy), are collected in MeshArray(meshes...). They are simply stored in a tuple. In cases when one wants to reduce the dimension of data, one can use MeshProduct(meshes...) to group meshes into a 1D flattened array.  

## Fourier Transform with DLR
We provide a dedicated fourier transform between imaginary time and Matsubara frequency spaces through the Discrete Lehmann representation.

```julia
    build green_dlr from above green's function (with to_dlr method)
    build green_tau from green_dlr (with dlr_to_tau method)
    Compare with exact analytical expression (can do this by plot)
    showcase the usage of <<, pipe in fourier transform
```
- The green's function in DLR space saves the spectral density of the object, which once generated, can be used to reproduce the Green's function in imaginary time/ Matsubara frequency space at any point with great accuracy (The reason why our current api always ask user to first go to DLR space)
- Explain the behavior of << for fourier transform, and when one should use it.


