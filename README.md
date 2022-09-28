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
We first show how to generate a multidimensional function given by the following equation: $G_{\sigma_1,\sigma_2}(q,\omega_n) = (Replace with G0 formula)$. This is a bare green's function of particles that have two inner spin 1/2 degrees, $sigma_1$ and $sigma_2$.  
    
TODO: Add the code block 
```julia
    build mesh of q
    build mesh of \omega_n
    build MeshArray
    Initialize by set value of data
    showcase intitialization with <<G0, G0 is a initialization function
```
- Explain how the inner spin states are represented in our structure
- Explain the q and \omega_n meshes forms a tuple 
- Explain the << operator

## Mesh Product
Now we would like to work with a more complicated example, where the momentum of $G_{\sigma_1,\sigma_2}(q,\omega_n)$ is two-dimensional. In this case we have to first build a 2D mesh as the direct product of two 1D meshes:

```julia
    build mesh of q1
    build mesh of q2
    build meshproduct of q1 and q2
    build mesh of \omega_n
    build MeshArray
    Initialize by set value of data
    showcase intitialization with <<G0, G0 is a initialization function
```
- Explain how the linear index mapping works in meshproduct

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


