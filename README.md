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

### Example 1: Green's function of a single level

We first show how to use `MeshArray` to present the Green's function of a single level quantum system filled with spinless fermionic particles. We assume that the system could exchange particles and energy with the enviroment and it is in an equilibrium state as a grand canonical ensemble. The single-particle Green's function is then has a simple form in Matsubara-frequency representation:  `G(ωₙ) = 1/(iωₙ - E)` where `E` is the level energy. We show how to generate and manipulate this Green's function.
     
```julia
julia> using GreenFunc, CompositeGrids

julia>  β = 100.0; E = 1.0 # inverse temperature and the level energy

julia> ωₙ_mesh = MeshGrids.ImFreq(100.0, FERMION; Euv = 10E) # UV energy cutoff is 10 times larger than the level energy
Matsubara frequency grid with 60 points, inverse temperature = 100.0, UV Energy scale = 10.0, fermionic = true: [-762, -374, -293, ..., 255, 394, 716]

julia> Gn =  MeshArray(ωₙ_mesh; dtype=ComplexF64); # Green's function defined on the ωₙ_mesh

julia> for (n, ωₙ) in enumerate(Gn.mesh[1])
           Gn[n] = 1/(ωₙ*im - E)
       end
```

- To define a Green's function, one first needs to define a set of meshes for the spacetime and other innerstates. This example demonstrate how to define a mesh in the temporal direction. The package provides three types of temporal meshes: imaginary time, matsubara frequnecy and DLR grids which provides a generic compressed representation for Green's functions (We will show how to use DLR later). They can be created with the `ImTime`, `ImFreq` and `DLRFreq` methods correspondingly. User needs to specify the inverse temperature, whether the particle is fermion or boson (using the constant `FERMION` or `BOSON`), and optionally the ultraviolet cutoff energy `Euv` and a set of predefined grid points `grid`. By default, a set of non-uniform grid points optimized for the given inverse temperature and the cutoff energy will be created.

- Once the meshes are created, one can define a `MeshArray` on them to represent the Green's function `Gn`. The constructor of `MeshArray` takes a set of meshes and initialze a multi-dimensional array. Each mesh corresponds to one dimension of the array. The data type of the `MeshArray` is specified by the optional keyword argument `dtype`, which is set to `Float64` by default. You can access the meshes (stored as a tuple) with `Gn.mesh`, and the array data with `Gn.data`.

- To initialize the `Gn.data`, one can use iterator interface of the meshes and the `MeshArray`. 

### Example 2: Green's function of a free electron gas

- [ ] Example use CompositeGrid to create the k-mesh with two inner states (spin up and down)

We first show how to generate a multidimensional function given by the following equation: $G_{ \sigma_1, \sigma_2}(q,\omega_n) = \frac{1}{i \omega_n - \epsilon}$. This is a bare green's function of particles that have two inner spin 1/2 degrees, $\sigma_1$ and $\sigma_2$, and a dispersion $\epsilon = q^2-1$.
     
```julia
julia> using GreenFunc, CompositeGrids

julia> kmesh = SimpleGrid.Uniform{Float64}([0.0, 10.0], 50); # initialze a uniform 

julia> mesh2 = MeshGrids.ImFreq(100.0, FERMION)
Matsubara frequency grid with 60 points, inverse temperature = 100.0, UV Energy scale = 10.0, fermionic = true: [-762, -374, -293, ..., 255, 394, 716]

julia> g_freq =  MeshArray(1:2, 1:2, mesh1, mesh2; dtype=ComplexF64) # Green's function with 2x2 innerstates
Meshed array with dims = (2, 2, 50, 60) and total length = 12000
- Mesh: Tuple{UnitRange{Int64}, UnitRange{Int64}, CompositeGrids.SimpleG.Uniform{Float64}, ImFreq{Float64, CompositeGrids.SimpleG.Arbitrary{Int64}}}

julia> for ind in eachindex(g_freq)
           q = g_freq.mesh[3][ind[3]]
           ω_n = g_freq.mesh[4][ind[4]]
           g_freq[ind] = 1/(im*ω_n - (q^2-1))
       end
```
- Momentum grids are handled by CompositeGrids package. Here SimpleGrid.Uniform{T}([mink,maxk],N) generates a linearly spaced momentum grid with N points in the interval (mink,maxk). For other specialized grids see https://github.com/numericalEFT/CompositeGrids.jl.
- GreenFunc package provides three specialized grids: imaginary time, matsubara frequnecy and DLR grids. They can be generated with ImTime, ImFreq and DLRFreq(beta,statistics) methods correspondingly. Here beta is the inverse temperature, and statistics defines whether the particle is fermion or boson. 
- By default the data type of MeshArray is Float64, and data is set to nothing if not provided. One can reset them with 'dtype' and 'data' keywards arguments. 
- Meshes for all variables, including the discrete ones (spin) and continuous ones (momentum and frequnecy), are collected in MeshArray(meshes...). They are simply stored in a tuple. In cases when one wants to reduce the dimension of data, one can use MeshProduct(meshes...) to group meshes into a 1D flattened array.  

### Example 3: Green's function of a Hubbard lattice

- [ ] Give an example use the uniform BZmesh to create the k-mesh

- [ ] Maybe discuss the MeshProduct here?

## Fourier Transform with DLR
We provide a dedicated fourier transform between imaginary time and Matsubara frequency spaces through the Discrete Lehmann representation.

```julia
#Transform the green function to DLR space and imaginary time space from Matsubara frequency spaces
julia > g_dlr = to_dlr(g)
julia > g_tau = dlr_to_imtime(g_dlr)

# Transform to Matsubara frequency space from DLR spectral density with "<<"
julia > g_freq_1 = MeshArray(innermesh..., mesh1, mesh2, dtype=ComplexF64)
julia > g_freq_1 << g_dlr

# Transform to imaginary time space from DLR spectral density with "<<"
julia > mesh3 = MeshGrids.ImTime(100.0, isFermi)
julia > g_tau_1 =  MeshArray(innermesh..., mesh1, mesh3, dtype=ComplexF64)
julia > g_tau_1 << g_dlr

# Transform to DLR space from imaginary time space with "<<"
julia > mesh4 = MeshGrids.DLRFreq(100,isFermi)
julia > g_dlr_1 =  MeshArray(innermesh..., mesh1, mesh4, dtype=ComplexF64)
julia > g_dlr_1 << g_tau

# julia > g_tau_0 =  MeshArray(innermesh..., mesh1, mesh3, dtype=ComplexF64)
# julia> for ind in eachindex(g_tau_0)
#            q = g.mesh[3][ind[3]]
#            τ = g.mesh[4][ind[4]]
#            g[ind] = -exp(-τ*(q^2-1))*(1-1/(exp(100*(q^2-1))+1))
#        end

# build green_dlr from above green's function (with to_dlr method)
# build green_tau from green_dlr (with dlr_to_tau method)
# Compare with exact analytical expression (can do this by plot)
# showcase the usage of <<, pipe in fourier transform
```
- The green's function in DLR space saves the spectral density of the object, which once generated, can be used to reproduce the Green's function in imaginary time/ Matsubara frequency space at any point with great accuracy (The reason why our current api always ask user to first go to DLR space)
- Explain the behavior of << for fourier transform, and when one should use it.


