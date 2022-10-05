# GreenFunc



[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://numericalEFT.github.io/GreenFunc.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://numericalEFT.github.io/GreenFunc.jl/dev)
[![Build Status](https://github.com/numericalEFT/GreenFunc.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/numericalEFT/GreenFunc.jl/actions/workflows/CI.yml?query=branch%3Amaster)
[![Coverage](https://codecov.io/gh/numericalEFT/GreenFunc.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/numericalEFT/GreenFunc.jl)

GreenFunc.jl differentiable numerical framework to manipulate multidimensional Green's functions.

## Features
 - `MeshArray` type provides an array defined on meshes. It provides a generic data structure for Green's functions.
 - Structured (non-)uniform temporal meshes for imaginary time/Matsubara frequency domains.
 - Compression and Fourier transform based on the discrete Lehmann representation.
 - Interface to the TRIQS library.
 
## Installation
This package has been registered. So, simply type `import Pkg; Pkg.add("GreenFunc")` in the Julia REPL to install.

## Basic Usage

### Example 1: Green's function of a single level

We first show how to use `MeshArray` to present the Green's function of a single level quantum system filled with spinless fermionic particles. We assume that the system could exchange particles and energy with the enviroment and it is in an equilibrium state as a grand canonical ensemble. The single-particle Green's function is then has a simple form in Matsubara-frequency representation:  $G(ωₙ) = \frac{1}{(iωₙ - E)}$ where `E` is the level energy. We show how to generate and manipulate this Green's function.
     
```julia
    using GreenFunc, CompositeGrids

    β = 100.0; E = 1.0 # inverse temperature and the level energy
    ωₙ_mesh = MeshGrids.ImFreq(100.0, FERMION; Euv = 100E) # UV energy cutoff is 100 times larger than the level energy
    Gn =  MeshArray(ωₙ_mesh; dtype=ComplexF64); # Green's function defined on the ωₙ_mesh

    for (n, ωₙ) in enumerate(Gn.mesh[1])
        Gn[n] = 1/(ωₙ*im - E)
    end
```

- To define a Green's function, one first needs to define a set of meshes for the spacetime and other innerstates. This example demonstrate how to define a mesh in the temporal direction. The package provides three types of temporal meshes: imaginary time, matsubara frequnecy and DLR grids which provides a generic compressed representation for Green's functions (We will show how to use DLR later). They can be created with the `ImTime`, `ImFreq` and `DLRFreq` methods correspondingly. User needs to specify the inverse temperature, whether the particle is fermion or boson (using the constant `FERMION` or `BOSON`), and optionally the ultraviolet cutoff energy `Euv` and a set of predefined grid points `grid`. By default, a set of non-uniform grid points optimized for the given inverse temperature and the cutoff energy will be created.

- Once the meshes are created, one can define a `MeshArray` on them to represent the Green's function `Gn`. The constructor of `MeshArray` takes a set of meshes and initialze a multi-dimensional array. Each mesh corresponds to one dimension of the array. The data type of the `MeshArray` is specified by the optional keyword argument `dtype`, which is set to `Float64` by default. You can access the meshes (stored as a tuple) with `Gn.mesh`, and the array data with `Gn.data`.

- By default `Gn.data` is left undefined if not specified by user. To initialize it, one can either use the  optional keyword argument `data` in constructor, or use iterator interface of the meshes and the `MeshArray`. 

### Example 2: Green's function of a free electron gas

Now let us show how to create a Green's function of a free electron gas. Compared with the single level spinless fermionic particles, free electrons have a spin-1/2 innerstates, and kinetic energy $ϵ_q = q^2-E$ (we use the unit where $m_e = 1/2$). The Green's function in Matsubara-frequnecy space is given by the following equation: $G_n = G_{\sigma_1, \sigma_2}(q,\omega_n) = \frac{1}{(i \omega_n - \epsilon_q)}$, where $\sigma_i$ denotes the spins of two electrons in the propagator. We inherit the Matsubara-frequency grid from the first example. We show how to use `CompositeGrids`package to generate momentum grids, and how multiple innerstates and meshes are treated by `MeshArray`.
```julia
    using GreenFunc, CompositeGrids
    β = 100.0; E = 1.0 # inverse temperature and the level energy
    ωₙ_mesh = MeshGrids.ImFreq(100.0, FERMION; Euv = 100E) # UV energy cutoff is 100 times larger than the level energy
    kmesh = SimpleGrid.Uniform{Float64}([0.0, 10.0], 50); # initialze an uniform momentum grid
    G_n =  MeshArray(1:2, 1:2, kmesh, ωₙ_mesh; dtype=ComplexF64); # Green's function of free electron gas with 2x2 innerstates

    for ind in eachindex(G_n)
        q = G_n.mesh[3][ind[3]]
        ω_n = G_n.mesh[4][ind[4]]
        G_n[ind] = 1/(im*ω_n - (q^2-E))
    end
```
- One can generate various types of grids with `CompositeGrids` package. The `SimpleGrid` module here provides several basic grids, such as uniform grids and logarithmically densed grids. The` Uniform` method here generates a 1D linearly spaced grid. User has to specify the number of grid points `N` and the boundary points `[min,max]`. One can also combine arbitrary numbers of `SimpleGrid` subgrids with a user specified pattern defined by a `panel grid`. These more advanced grids optimized for different purposes can be found on https://github.com/numericalEFT/CompositeGrids.jl.

- The constructor of `MeshArray` can take any iterable objects as one of its meshes. Therefore for discrete innerstates such as spins, one can simply use an integer tuple.

### Example 3: Green's function of a Hubbard lattice

Now we show how to generate a multidimensional function with
tight-binding dispersion:
$G(q, \omega_n) = \frac{1}{i\omega_n - \epsilon_q}$ with 
$\epsilon_q = -2t(cos(q_x)+cos(q_y))$.
The momentum is defined on the first Brillouin zone captured by a 2D k-mesh.

```julia
    using GreenFunc
    using GreenFunc: BrillouinZoneMeshes

    DIM, nk = 2, 8
    latvec = [1.0 0.0; 0.0 1.0] .* 2π
    bzmesh = BrillouinZoneMeshes.BaseMesh.UniformMesh{DIM, nk}([0.0, 0.0], latvec)
    ωₙmesh = MeshGrids.ImFreq(10.0, FERMION)
    g_freq =  MeshArray(bzmesh, ωₙmesh; dtype=ComplexF64)

    t = 1.0
    for ind in eachindex(g_freq)
        q = g_freq.mesh[1][ind[1]]
        ωₙ = g_freq.mesh[2][ind[2]]
        g_freq[ind] = 1/(ωₙ*im - (-2*t*sum(cos.(q))))
    end
    g_tau = g_freq |> to_dlr |> to_imtime #fourier transform to (k, tau) domain
```

- For lattice systems with multidimenstional Brillouin zone, the momentum grids can be better handled with the BrillouinZoneMeshes package. Here a `UniformMesh{DIM,N}(origin, latvec)` generates a linearly spaced momentum mesh on the first Brillouin zone defined by origin and lattice vectors given. For more detail see https://github.com/numericalEFT/BrillouinZoneMeshes.jl.


### Example 4:  Fourier Transform of Green's function with DLR
In this example we show one of the key features of our package: a dedicated fourier transform between imaginary time and Matsubara frequency spaces through the Discrete Lehmann representation (DLR). The DLR allows compact and efficient storage of any functions that has Lehmann representation. It saves just enough information in spectral density on optimized real frequency points to represent the function up to user specified accuracy `ϵ`. See https://github.com/numericalEFT/Lehmann.jl for more details.

Let us first convert our free electron Green's function from Matsubara frequency space `G_n` to DLR space `G_dlr`, and then convert `G_dlr` back to `G_n2` that is identical to `G_n`.
```julia
    using GreenFunc, CompositeGrids

    β = 100.0; E = 1.0 # inverse temperature and the level energy
    ωₙ_mesh = MeshGrids.ImFreq(100.0, FERMION; Euv = 100E) # UV energy cutoff is 100 times larger than the level energy
    kmesh = SimpleGrid.Uniform{Float64}([0.0, 10.0], 50); # initialze an uniform momentum grid
    G_n =  MeshArray(1:2, 1:2, kmesh, ωₙ_mesh; dtype=ComplexF64); # Green's function of free electron gas with 2x2 innerstates

    for ind in eachindex(G_n)
        q = G_n.mesh[3][ind[3]]
        ω_n = G_n.mesh[4][ind[4]]
        G_n[ind] = 1/(im*ω_n - (q^2-E))
    end

    G_dlr = to_dlr(G_n) # convert G_n to DLR space 
    G_n2 = dlr_to_imfreq(G_dlr) # convert G_dlr back to Matsubara frequency space
    @assert G_n2 ≈ G_n

```
For the second step, let us construct the same Green's function in imaginary time space given by $G_{\tau} = -e^{-\tau \epsilon_q}/(1+e^{-\beta \epsilon_q})$, and reproduce it by transforming `G_dlr`.
```julia
    τ_mesh = MeshGrids.ImTime(100.0, FERMION; Euv = 100E) 
    G_τ =  MeshArray(1:2, 1:2, kmesh, τ_mesh; dtype=Float64)

    for ind in eachindex(G_τ)
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

    G_τ2 = dlr_to_imtime(G_dlr,τ_mesh) #Transform G_dlr to imaginary time space.
    @assert G_τ ≈ G_τ2

    G_τ3 = G_n |> to_dlr |> dlr_to_imtime; #Use pipe to do a two-step transform
    @assert G_τ ≈ G_τ3
```

- For any Green's function that has eaxctly one temperal grid (`ImTime`, `ImFreq` and `DLRFreq`) in meshes, we provide a set of operations (`to_dlr`, `dlr_to_imfreq` and `dlr_to_imtime`) to bridge the DLR space with imaginary-time and Matsubara-frequency space. By default, all of them find the dimension of the temperal grid within meshes, and perform the transformation with respect to it. Otherwise, one can always specify the dimension with the optional keyword argumen `dim`. 

- The `to_dlr` method calculate the DLR spectral density of a Green's function in imaginary-time or Matsubara-frequency space, and save it in a new `MeshArray` object. This operation requires a pre-defined `DLRFreq` grid that has the same inverse temperature and statistics (fermion or boson) as the `ImFreq` or `ImTime` grids in the original Green's function. By default, the DLR grid is automatically generated with accuracy `ϵ=1e-12`. User can also specify the accuracy with keywards argument `rtol`, or provide the DLR grid.

- Once a spectral density `G_dlr` in DLR space is obtained, one can use `dlr_to_imfreq` or `dlr_to_imtime` methods to reconstruct the Green's function in corresponding space. Since DLR saves information for the entire function, the target `ImFreq` or `ImTime` grid can be any grids specified by user. By default it is set to be the optimized grid from `DLRFreq` mesh in `G_dlr`. 

- Since the spectral density `G_dlr` can be reused whenever user wants to change the grid points of Green's function (normally through interpolation that lost more accuracy then DLR transfrom), we encourage user to always keep the `G_dlr` objects. User can use piping operator |> as shown to do fourier transform directly between `ImFreq` and `ImTime` in one line, although it will throw away the spectral density.

##  Interface with TRIQS

TRIQS (Toolbox for Research on Interacting Quantum Systems) is a scientific project providing a set of C++ and Python libraries for the study of interacting quantum systems. We provides direct interface to convert TRIQS objects, such as temporal meshes, Brillioun zone meshes and multi-dimensional Green's functions, to the equivalent objects in our package. This helps TRIQS user to make use of our package without worrying about the different interal data structures.

### Example 5: Load Triqs Temporal Mesh
First we show how to import an imaginary time mesh from TRIQS.
```julia
    using PythonCall, Greenfunc
    gf = pyimport("triqs.gf")
    np = pyimport("numpy")

    mt = gf.MeshImTime(beta=1.0, S="Fermion", n_max=3)
    mjt = from_triqs(mt)
    for (i, x) in enumerate([p for p in mt.values()])
        @assert mjt[i] ≈ pyconvert(Float64, x)
    end
    
```
- With the `PythonCall` package, one can import python packages with `pyimport` and directly exert python code in julia. Here we import the green's function module `triqs.gf`, and generate an uniform imaginary time mesh with `MeshImTime`. User has to specify the inverse temperature,  whether the particle is fermion or boson and the number of grid points.

- Once a TRIQS object is prepared, one can simply convert it to an equivalent object in our package with `from_triqs`. The object can be a mesh, a Green's function or a block Green's function. In this example the TRIQS imaginary time grid is converted to an identical `ImTime` grid.

- The objects generated by python code have to be converted to corresponding Julia types with `pyconvert`. Here all values in the TRIQS mesh are python float numbers. After they are converted to Julia type Float64, they can be directly compared with grid points in `ImTime` grid we generated. 
### Example 6: Load Triqs BrillouinZone

In this example we show how the Brillouin zone mesh from triqs is loaded
as a UniformMesh from BrillouinZoneMeshes package, and clarify the
convention we adopted to resolve the conflict between Julia and Python
data structure.

```julia
    using PythonCall, GreenFunc

    # construct triqs Brillouin zone mesh
    lat = pyimport("triqs.lattice")
    gf = pyimport("triqs.gf")
    BL = lat.BravaisLattice(units=((2, 0, 0), (1, sqrt(3), 0))) 
    BZ = lat.BrillouinZone(BL)
    nk = 4
    mk = gf.MeshBrillouinZone(BZ, nk)

    # load Triqs mesh and construct 
    mkj = from_triqs(mk)

    for p in mk
        pval = pyconvert(Array, p.value)
        # notice that Triqs always return a 3D point, even for 2D case(where z is always 0)
        # notice also that Julia index starts from 1 while Python from 0
        # points of the same linear index has the same value
        ilin = pyconvert(Int, p.linear_index) + 1
        @assert pval[1:2] ≈ mkj[ilin]
        # points with the same linear index corresponds to REVERSED cartesian index
        inds = pyconvert(Array, p.index)[1:2] .+ 1
        @assert pval[1:2] ≈ mkj[reverse(inds)...]
    end
```

- Due to the difference of cartesian index order between Julia and Python,
it is inconvenient to construct Brillouin zone mesh exactly same as in Triqs.
- We adopted the convention so that the grid point and linear index
is consistent with Triqs counterpart, while the order of cartesian index
and lattice vector reversed.
- Here's a table of 2D converted mesh v.s. the Triqs counterpart:

| Object          | Triqs             | GreenFunc.jl   |
| --------------- | ----------------- | -------------- |
| Linear index    | mk[i]=(x, y, 0)   | mkj[i]= (x, y) |
| Cartesian index | mk[i,j]=(x, y, 0) | mkj[j,i]=(x,y) |
| Lattice vector  | (a1, a2)          | (a2, a1)       |

### Example 7: Load Triqs Greens function of a Hubbard Lattice

A TRIQS Green's function is defined on a set of meshes of continuous variables, together with the discrete innerstates specified by the `target_shape`. The structure is immediately representable by `MeshArray`. In the following example, we show how to generate a TRIQS Green's function of a Hubbard lattice within Julia. The same as Example 3, the Green's function is given by $G(q, \omega_n) = \frac{1}{i\omega_n - \epsilon_q}$ with $\epsilon_q = -2t(cos(q_x)+cos(q_y))$. 

```julia
    using PythonCall, GreenFunc
    
    np = pyimport("numpy")
    lat = pyimport("triqs.lattice")
    gf = pyimport("triqs.gf")
    
    BL = lat.BravaisLattice(units=((2, 0, 0), (1, sqrt(3), 0))) # testing with a triangular lattice so that exchanged index makes a difference
    BZ = lat.BrillouinZone(BL)
    nk = 20
    mk = gf.MeshBrillouinZone(BZ, nk)
    miw = gf.MeshImFreq(beta=1.0, S="Fermion", n_max=100)
    mprod = gf.MeshProduct(mk, miw)

    G_w = gf.GfImFreq(mesh=miw, target_shape=[1, 1]) #G_w.data.shape will be [201, 1, 1]
    G_k_w = gf.GfImFreq(mesh=mprod, target_shape = [2, 3] ) #target_shape = [2, 3] --> innerstate = [3, 2]

    # Due to different cartesian index convention in Julia and Python, the data g_k_w[n, m, iw, ik] corresponds to G_k_w.data[ik-1, iw-1, m-1, n-1])

    t = 1.0
    for (ik, k) in enumerate(G_k_w.mesh[0])
        G_w << gf.inverse(gf.iOmega_n - 2 * t * (np.cos(k[0]) + np.cos(k[1])))
        G_k_w.data[ik-1, pyslice(0, nk^2), pyslice(0, G_k_w.target_shape[0]) , pyslice(0,G_k_w.target_shape[1])] = G_w.data[pyslice(0, nk^2), pyslice(0, G_w.target_shape[0]) , pyslice(0,G_w.target_shape[1])] #pyslice = :      
    end

    g_k_w = from_triqs(G_k_w)
    
    g_k_w2 = MeshArray(G_k_w) 
    @assert g_k_w2 ≈ g_k_w
    g_k_w2 << G_k_w
    @assert g_k_w2 ≈ g_k_w
    
```
- The `MeshProduct` from TRIQS is decomposed into separate meshes, and converted to corresponding meshes in our package to construct the `MeshArray`.
- The `target_shape` in TRIQS Green's function are converted to integer tuples that represents the discrete degrees of freedoms.
` Data slicing with `:` is not available in julia. One has to use pyslice instead.
- As explained in Example 6, the cartesian index order of data has to be inversed during the conversion.
- We support three different interfaces for conversion of TRIQS Green's function. One can construct a new MeshArray with `from_triqs` or `MeshArray` constructor. One can also load TRIQS Green's function into an existing `MeshArray` with the `<<` operator.

### Example 8: Load Triqs block Greens function

The block Greens function in TRIQS is a dictionary of different Green's functions. Below we show an example that converts it into a dictionary of `MeshArray` objects.

```julia
    using PythonCall, GreenFunc

    gf = pyimport("triqs.gf")
    np = pyimport("numpy")
    mt = gf.MeshImTime(beta=1.0, S="Fermion", n_max=3)
    lj = pyconvert(Int, @py len(mt))
    G_t = gf.GfImTime(mesh=mt, target_shape=[2,3]) #target_shape = [2, 3] --> innerstate = [3, 2]
    G_w = gf.GfImTime(mesh=mt, target_shape=[2,3]) #target_shape = [2, 3] --> innerstate = [3, 2]

    blockG = gf.BlockGf(name_list=["1", "2"], block_list=[G_t, G_w], make_copies=false)

    jblockG = from_triqs(blockG) 
    #The converted block Green's function is a dictionary of MeshArray corresponding to TRIQS block Green's function. The mapping between them is: jblockG["name"][i1, i2, t] = blockG["name"].data[t-1, i2-1, i1-1]

```

