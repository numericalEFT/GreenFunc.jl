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

julia> ωₙ_mesh = MeshGrids.ImFreq(100.0, FERMION; Euv = 100E) # UV energy cutoff is 100 times larger than the level energy
Matsubara frequency grid with 83 points, inverse temperature = 100.0, UV Energy scale = 100.0, fermionic = true: [-8083, -4302, -2816, ..., 3233, 4145, 8396]


julia> Gn =  MeshArray(ωₙ_mesh; dtype=ComplexF64); # Green's function defined on the ωₙ_mesh

julia> for (n, ωₙ) in enumerate(Gn.mesh[1])
           Gn[n] = 1/(ωₙ*im - E)
       end
```

- To define a Green's function, one first needs to define a set of meshes for the spacetime and other innerstates. This example demonstrate how to define a mesh in the temporal direction. The package provides three types of temporal meshes: imaginary time, matsubara frequnecy and DLR grids which provides a generic compressed representation for Green's functions (We will show how to use DLR later). They can be created with the `ImTime`, `ImFreq` and `DLRFreq` methods correspondingly. User needs to specify the inverse temperature, whether the particle is fermion or boson (using the constant `FERMION` or `BOSON`), and optionally the ultraviolet cutoff energy `Euv` and a set of predefined grid points `grid`. By default, a set of non-uniform grid points optimized for the given inverse temperature and the cutoff energy will be created.

- Once the meshes are created, one can define a `MeshArray` on them to represent the Green's function `Gn`. The constructor of `MeshArray` takes a set of meshes and initialze a multi-dimensional array. Each mesh corresponds to one dimension of the array. The data type of the `MeshArray` is specified by the optional keyword argument `dtype`, which is set to `Float64` by default. You can access the meshes (stored as a tuple) with `Gn.mesh`, and the array data with `Gn.data`.

- By default `Gn.data` is left undefined if not specified by user. To initialize it, one can either use the  optional keyword argument `data` in constructor, or use iterator interface of the meshes and the `MeshArray`. 

### Example 2: Green's function of a free electron gas

- [ ] Example use CompositeGrid to create the k-mesh with two inner states (spin up and down)

Now let us show how to create a Green's function of a free electron gas. Compared with the single level spinless fermionic particles, free electrons have a spin-1/2 innerstates, and kinetic energy `ϵ_q = q^2-E` (we use the unit where `m_e = 1/2`). The Green's function in Matsubara-frequnecy space is given by the following equation: `G_n = G_{ σ_1, σ_2}(q,ω_n) = 1/(i ω_n - ϵ_q)`, where `σ_i` denotes the spins of two electrons in the propagator. We inherit the Matsubara-frequency grid from the first example. We show how to use `CompositeGrids`package to generate momentum grids, and how multiple innerstates and meshes are treated by `MeshArray`.
```julia
julia> kmesh = SimpleGrid.Uniform{Float64}([0.0, 10.0], 50); # initialze an uniform momentum grid

julia> G_n =  MeshArray(1:2, 1:2, kmesh, ωₙ_mesh; dtype=ComplexF64); # Green's function of free electron gas with 2x2 innerstates

julia> for ind in eachindex(G_n)
           q = G_n.mesh[3][ind[3]]
           ω_n = G_n.mesh[4][ind[4]]
           G_n[ind] = 1/(im*ω_n - (q^2-E))
       end
```
- One can generate various types of grids with `CompositeGrids` package. The `SimpleGrid` module here provides several basic grids, such as uniform grids and logarithmically densed grids. The` Uniform` method here generates a 1D linearly spaced grid. User has to specify the number of grid points `N` and the boundary points `[min,max]`. One can also combine arbitrary numbers of `SimpleGrid` subgrids with a user specified pattern defined by a `panel grid`. These mored advanced grids optimized for different purposes can be found on https://github.com/numericalEFT/CompositeGrids.jl.

- The constructor of `MeshArray` can take any iterable objects as one of its meshes. Therefore for discrete innerstates such as spins, one can simply use an integer tuple.

### Example 3: Green's function of a Hubbard lattice

- [ ] Give an example use the uniform BZmesh to create the k-mesh

Now we show how to generate a multidimensional function with
tight-binding dispersion:
$G_(k, \omega_n) = \frac{1}{i\omega_n - \epsilon(k)}$ with 
$\epsilon(k) = -2t(cos(k_x)+cos(k_y))$.
The momentum is defined on the first Brillouin zone captured by a 2D k-mesh.

```julia
using GreenFunc
using GreenFunc: BrillouinZoneMeshes

DIM = 2
latvec = [1.0 0.0; 0.0 1.0] .* 2π
nk = 8
bzmesh = BrillouinZoneMeshes.BaseMesh.UniformMesh{DIM, nk}([0.0, 0.0], latvec)
wmesh = MeshGrids.ImFreq(10.0, FERMION)
g_freq =  MeshArray(bzmesh, wmesh; dtype=ComplexF64)

t = 1.0
for ind in eachindex(g_freq)
    q = g_freq.mesh[1][ind[1]]
    ω_n = g_freq.mesh[2][ind[2]]
    println(q, ω_n)
    g_freq[ind] = 1/(im*ω_n - (-2*t*sum(cos.(q))))
end
```

- Momentum is handled by BrillouinZoneMeshes package. Here a UniformMesh{DIM,N}(origin, latvec) generates a linearly spaced momentum mesh on the first Brillouin zone defined by origin and lattice vectors given. For more detail see https://github.com/numericalEFT/BrillouinZoneMeshes.jl.

- [ ] Maybe discuss the MeshProduct here?

### Example 4:  Fourier Transform of Green's function with DLR
In this example we show one of the key features of our package: a dedicated fourier transform between imaginary time and Matsubara frequency spaces through the Discrete Lehmann representation (DLR). The DLR allows compact and efficient storage of any functions that has Lehmann representation. It saves just enough information in spectral density on optimized real frequency points to represent the function up to user specified accuracy `ϵ`. See https://github.com/numericalEFT/Lehmann.jl for more details.

Let us first convert our free electron Green's function from Matsubara frequency space `G_n` to DLR space `G_dlr`, and then convert `G_dlr` back to `G_n2` that is identical to `G_n`.
```julia
julia> G_dlr = to_dlr(G_n) # convert G_n to DLR space 
Meshed array with dims = (2, 2, 50, 60) and total length = 12000
- Mesh: Tuple{UnitRange{Int64}, UnitRange{Int64}, CompositeGrids.SimpleG.Uniform{Float64}, DLRFreq{Float64}} 

julia> G_n2 = dlr_to_imfreq(G_dlr) # convert G_dlr back to Matsubara frequency space
Meshed array with dims = (2, 2, 50, 60) and total length = 12000
- Mesh: Tuple{UnitRange{Int64}, UnitRange{Int64}, CompositeGrids.SimpleG.Uniform{Float64}, ImFreq{Float64, CompositeGrids.SimpleG.Arbitrary{Int64}}} 

julia> G_n2 ≈ G_n
true
```
For the second step, let us construct the same Green's function in imaginary time space given by `G_τ = -e^{-τϵ_q}/(1+e^{-βϵ_q})`, and reproduce it by transforming `G_dlr`.
```julia
julia> τ_mesh = MeshGrids.ImTime(100.0, FERMION; Euv = 100E) 
Imaginary Time grid with 83 points, inverse temperature = 100.0, UV Energy scale = 100.0, fermionic = true: [2.6136191e-5, 0.0020572557, 0.0068079994, ..., 99.993192, 99.997943, 99.999974]

julia> G_τ =  MeshArray(1:2, 1:2, kmesh, τ_mesh; dtype=Float64);

julia> for ind in eachindex(G_τ)
           q = G_τ.mesh[3][ind[3]]
           τ = G_τ.mesh[4][ind[4]]
           ϵ_q = q^2-E
           x = β*ϵ_q/2.0
           y = 2*τ/β - 1.0
           if x>0.0 #Make sure the index is always negative to avoid overflow issue from exponential function
               G_τ[ind] =-exp(-x*(y+1))/(1+exp(-2*x))
           else
               G_τ[ind] =-exp(-x*(y-1))/(1+exp(2*x))
           end
       end

julia> G_τ2 = dlr_to_imtime(G_dlr,τ_mesh) #Transform G_dlr to imaginary time space.

julia> G_τ ≈ G_τ2
true

julia> G_τ3 = G_n |> to_dlr |> dlr_to_imtime; #Use pipe to do a two-step transform

julia> G_τ ≈ G_τ3
true

```

- For any Green's function that has eaxctly one temperal grid (`ImTime`, `ImFreq` and `DLRFreq`) in meshes, we provide a set of operations (`to_dlr`, `dlr_to_imfreq` and `dlr_to_imtime`) to bridge the DLR space with imaginary-time and Matsubara-frequency space. By default, all of them find the dimension of the temperal grid within meshes, and perform the transformation with respect to it. Otherwise, one can always specify the dimension with the optional keyword argumen `dim`. 

- The `to_dlr` method calculate the DLR spectral density of a Green's function in imaginary-time or Matsubara-frequency space, and save it in a new `MeshArray` object. This operation requires a pre-defined `DLRFreq` grid that has the same inverse temperature and statistics (fermion or boson) as the `ImFreq` or `ImTime` grids in the original Green's function. By default, the DLR grid is automatically generated with accuracy `ϵ=1e-12`. User can also specify the accuracy with keywards argument `rtol`, or provide the DLR grid.

- Once a spectral density `G_dlr` in DLR space is obtained, one can use `dlr_to_imfreq` or `dlr_to_imtime` methods to reconstruct the Green's function in corresponding space. Since DLR saves information for the entire function, the target `ImFreq` or `ImTime` grid can be any grids specified by user. By default it is set to be the optimized grid from `DLRFreq` mesh in `G_dlr`. 

- Since the spectral density `G_dlr` can be reused whenever user wants to change the grid points of Green's function (normally through interpolation that lost more accuracy then DLR transfrom), we encourage user to always keep the `G_dlr` objects. User can use piping operator |> as shown to do fourier transform directly between `ImFreq` and `ImTime` in one line, although it will throw away the spectral density.

##  Triqs

### Example 5: Load Triqs Temporal Mesh
Another 

### Example 6: Load Triqs BrillouinZone


### Example 7: Load Triqs Green;s function
