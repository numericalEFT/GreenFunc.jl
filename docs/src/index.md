```@meta
CurrentModule = GreenFunc
```

# GreenFunc



[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://numericalEFT.github.io/GreenFunc.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://numericalEFT.github.io/GreenFunc.jl/dev)
[![Build Status](https://github.com/numericalEFT/GreenFunc.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/numericalEFT/GreenFunc.jl/actions/workflows/CI.yml?query=branch%3Amaster)
[![Coverage](https://codecov.io/gh/numericalEFT/GreenFunc.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/numericalEFT/GreenFunc.jl)

GreenFunc.jl is a differentiable numerical framework to manipulate multidimensional Green's functions.

## Features
 - `MeshArray` type as an array defined on meshes, which provides a generic data structure for Green's functions, vertex functions or any other correlation/response functions.
 - Structured (non-)uniform Brillouin Zone meshes powered by the package [`BrillouinZoneMeshes.jl`](https://github.com/numericalEFT/BrillouinZoneMeshes.jl).
 - Structured (non-)uniform temporal meshes for (imaginary-)time or (Matsubara-)frequency domains powered by the pacakge [`CompositeGrids.jl`](https://github.com/numericalEFT/CompositeGrids.jl).
 - Compat representation based on the Discrete Lehmann representation (DLR) powered by the package [`Lehmann.jl`](https://github.com/numericalEFT/Lehmann.jl).
 - Accurate and fast Fourier transform.
 - Interface to the [`TRIQS`](https://triqs.github.io/) library.
 
## Installation
This package has been registered. So, simply type `import Pkg; Pkg.add("GreenFunc")` in the Julia REPL to install.

## Basic Usage

### Example 1: Green's function of a single level

We first show how to use `MeshArray` to present Green's function of a single-level quantum system filled with spinless fermionic particles. We assume that the system could exchange particles and energy with the environment so that it's equilibrium state is a grand canonical ensemble. The single-particle Green's function then has a simple form in Matsubara-frequency representation:  $G(ωₙ) = \frac{1}{(iωₙ - E)}$ where $E$ is the level energy. We show how to generate and manipulate this Green's function.

```julia
using GreenFunc

# inverse temperature and the level energy
β = 100.0; E = 1.0
# UV energy cutoff is 100 times larger than the level energy
ωₙ_mesh = MeshGrids.ImFreq(100.0, FERMION; Euv = 100E)
# Green's function defined on the ωₙ_mesh
Gn = MeshArray(ωₙ_mesh; dtype=ComplexF64)

for (n, ωₙ) in enumerate(Gn.mesh[1])
    Gn[n] = 1/(ωₙ*im - E)
end
```

- Green's function describes correlations between two or more spacetime events. The spacetime continuum needs to be discretized into spatial and temporal meshes. This example demonstrates how to define a one-body Green's function on a temporal mesh. The package provides three types of temporal meshes: imaginary-time grid, Matsubara-frequency grid, and DLR grid. The latter provides a generic compressed representation for Green's functions (We will show how to use DLR later).  Correspondingly, They can be created with the `ImTime`, `ImFreq`, and `DLRFreq` methods. The user needs to specify the inverse temperature, whether the particle is fermion or boson (using the constant `FERMION` or `BOSON`). Internally, a set of non-uniform grid points optimized for the given inverse temperature and the cutoff energy will be created with the given parameters.

- Once the meshes are created, one can define a `MeshArray` on them to represent the Green's function `Gn`. The constructor of `MeshArray` takes a set of meshes and initializes a multi-dimensional array. Each mesh corresponds to one dimension of the array. The data type of the `MeshArray` is specified by the optional keyword argument `dtype`, which is set to `Float64` by default. You can access the meshes (stored as a tuple) with `Gn.mesh` and the array data with `Gn.data`.

- By default, `Gn.data` is left undefined if not specified by the user. To initialize it, one can either use the optional keyword argument `data` in the constructor or use the iterator interface of the meshes and the `MeshArray`. 

### Example 2: Green's function of a free electron gas

Now let us show how to create a Green's function of a free electron gas. Unlike the spinless fermionic particle, the electron is a spin-1/2 particle so that it has two inner states. In free space, it has a kinetic energy $ϵ_q = q^2-E$ (we use the unit where $m_e = 1/2$). The Green's function in Matsubara-frequency space is then given by the following equation: $G_n = G_{\sigma_1, \sigma_2}(q,\omega_n) = \frac{1}{i \omega_n - \epsilon_q}$, where $\sigma_i$ denotes the spins of the incoming and the outgoing electron in the propagator. We inherit the Matsubara-frequency grid from the first example. We show how to use the `CompositeGrids` package to generate momentum grids and how to treat the multiple inner states and the meshes with `MeshArray`.
```julia
using GreenFunc, CompositeGrids

# inverse temperature and the level energy
β = 100.0; E = 1.0
# UV energy cutoff is 100 times larger than the level energy
ωₙ_mesh = MeshGrids.ImFreq(100.0, FERMION; Euv = 100E)
# initialze an uniform momentum grid
kmesh = SimpleGrid.Uniform{Float64}([0.0, 10.0], 50)
# Green's function of free electron gas with 2x2 innerstates
G_n =  MeshArray(1:2, 1:2, kmesh, ωₙ_mesh; dtype=ComplexF64)

for ind in eachindex(G_n)
    q = G_n.mesh[3][ind[3]]
    ω_n = G_n.mesh[4][ind[4]]
    G_n[ind] = 1/(ω_n*im - (q^2-E))
end
```
- One can generate various types of grids with the `CompositeGrids` package. The `SimpleGrid` module here provides several basic grids, such as uniform grids and logarithmically dense grids. The` Uniform` method here generates a 1D linearly spaced grid. The user has to specify the number of grid points `N` and the boundary points `[min, max]`. One can also combine arbitrary numbers of `SimpleGrid` subgrids with a user-specified pattern defined by a `panel grid`. These more advanced grids optimized for different purposes can be found in this [link](https://github.com/numericalEFT/CompositeGrids.jl).

- The constructor of `MeshArray` can take any iterable objects as one of its meshes. Therefore for discrete inner states such as spins, one can simply use a `1:2`, which is a `UnitRange{Int64}` object.

### Example 3: Green's function of a Hubbard lattice

Now we show how to generate a multi-dimensional Green's function on a Brillouin Zone meshe. We calculate the Green's function of a free spinless Fermi gas on a square lattice. It has a tight-binding dispersion $\epsilon_q = -2t(\cos(q_x)+\cos(q_y))$, which gives
$G(q, \omega_n) = \frac{1}{i\omega_n - \epsilon_q}$.
The momentum is defined on the first Brillouin zone captured by a 2D k-mesh.

```julia
using GreenFunc
using GreenFunc: BrillouinZoneMeshes

DIM, nk = 2, 8
latvec = [1.0 0.0; 0.0 1.0] .* 2π
bzmesh = BrillouinZoneMeshes.BaseMesh.UniformMesh{DIM, nk}([0.0, 0.0], latvec)
ωₙmesh = ImFreq(10.0, FERMION)
g_freq =  MeshArray(bzmesh, ωₙmesh; dtype=ComplexF64)

t = 1.0
for ind in eachindex(g_freq)
    q = g_freq.mesh[1][ind[1]]
    ωₙ = g_freq.mesh[2][ind[2]]
    g_freq[ind] = 1/(ωₙ*im - (-2*t*sum(cos.(q))))
end
```

- For lattice systems with multi-dimensional Brillouin zone, the momentum grids internally generated with the [`BrillouinZoneMeshes.jl`](https://github.com/numericalEFT/BrillouinZoneMeshes.jl) package. Here a `UniformMesh{DIM,N}(origin, latvec)` generates a linearly spaced momentum mesh on the first Brillouin zone defined by origin and lattice vectors given. For more detail, see the [link](https://github.com/numericalEFT/BrillouinZoneMeshes.jl).


### Example 4:  Fourier Transform of Green's function with DLR
DLR provides a compact representation for one-body Green's functions. At a temperature $T$ and an accuracy level $\epsilon$, it represents a generic Green's function with only $\log (1/T) \log (1/\epsilon)$ basis functions labeled by a set of real frequency grid points. It enables fast Fourier transform and interpolation between the imaginary-time and the Matsubara-frequency representations with a cost $O(\log (1/T) \log (1/\epsilon))$. `GreenFunc.jl` provide DLR through the package [`Lehmann.jl`](https://github.com/numericalEFT/Lehmann.jl).

In the following example, we demonstrate how to perform DLR-based Fourier transform in `GreenFunc.jl` between the imaginary-time and the Matsubara-frequency domains back and forth through the DLR representation.
```julia
using GreenFunc, CompositeGrids

β = 100.0; E = 1.0 # inverse temperature and the level energy
ωₙ_mesh = ImFreq(100.0, FERMION; Euv = 100E) # UV energy cutoff is 100 times larger than the level energy
kmesh = SimpleGrid.Uniform{Float64}([0.0, 10.0], 50); # initialze an uniform momentum grid
G_n =  MeshArray(1:2, 1:2, kmesh, ωₙ_mesh; dtype=ComplexF64); # Green's function of free electron gas with 2x2 innerstates

for ind in eachindex(G_n)
    q = G_n.mesh[3][ind[3]]
    ω_n = G_n.mesh[4][ind[4]]
    G_n[ind] = 1/(im*ω_n - (q^2-E))
end

G_dlr = to_dlr(G_n) # convert G_n to DLR space
G_tau = to_imtime(G_dlr) # convert G_dlr to the imaginary-time domain

#alternative, you can use the pipe operator
G_tau = G_n |> to_dlr |> to_imtime #Fourier transform to (k, tau) domain

```
The imaginary-time Green's function after the Fourier transform shoud be consistent with the analytic solution $G_{\tau} = -e^{-\tau \epsilon_q}/(1+e^{-\beta \epsilon_q})$.

- For any Green's function that has at least one imaginary-temporal grid (`ImTime`, `ImFreq`, and `DLRFreq`) in meshes, we provide a set of operations (`to_dlr`, `to_imfreq` and `to_imtime`) to bridge the DLR space with imaginary-time and Matsubara-frequency space. By default, all these functions find the dimension of the imaginary-temporal mesh within Green's function meshes and perform the transformation with respect to it. Alternatively, one can specify the dimension with the optional keyword argument `dim`. Be careful that the original version of DLR is only guaranteed to work with one-body Green's function.

- Once a spectral density `G_dlr` in DLR space is obtained, one can use `to_imfreq` or `to_imtime` methods to reconstruct the Green's function in the corresponding space. By default, `to_imfreq` and `to_imtime` uses an optimized imaginary-time or Matsubara-frequency grid from the DLR. User can assign a target imaginary-time or Matsubara-frequency grid if necessary.   

- Combining `to_dlr`, `to_imfreq` and `to_imtime` allows both _interpolation_ as well as _Fourier transform_.

- Since the spectral density `G_dlr` can be reused whenever the user wants to change the grid points of Green's function (normally through interpolation that lost more accuracy than DLR transform), we encourage the user always to keep the `G_dlr` objects. If the intermediate DLR Green's function is not needed, the user can use piping operator `|>` as shown to do Fourier transform directly between `ImFreq` and `ImTime` in one line.


##  Interface with TRIQS

TRIQS (Toolbox for Research on Interacting Quantum Systems) is a scientific project providing a set of C++ and Python libraries for the study of interacting quantum systems. We provide a direct interface to convert TRIQS objects, such as the temporal meshes, the Brillouin zone meshes, and the  multi-dimensional (blocked) Green's functions, to the equivalent objects in our package. It would help TRIQS users to make use of our package without worrying about the different internal data structures.

The interface is provided by an independent package [`NEFTInterface.jl`](https://github.com/numericalEFT/NEFTInterface.jl). We provide several examples of interfacing TRIQS and `GreenFunc.jl` in the [`NEFTInterface.jl` README](https://github.com/numericalEFT/NEFTInterface.jl).
