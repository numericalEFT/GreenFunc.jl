var documenterSearchIndex = {"docs":
[{"location":"lib/mesharrays/#Array-with-MeshGrids","page":"MeshArrays","title":"Array with MeshGrids","text":"","category":"section"},{"location":"lib/mesharrays/","page":"MeshArrays","title":"MeshArrays","text":"Modules = [GreenFunc.MeshArrays]","category":"page"},{"location":"lib/mesharrays/#GreenFunc.MeshArrays.MeshArray","page":"MeshArrays","title":"GreenFunc.MeshArrays.MeshArray","text":"struct MeshArray{T,N,MT} <: AbstractMeshArray{T,N}\n\nMulti-dimensional array that is defined on a mesh.  The mesh is a tuple of meshgrid objects.  The mesh is stored in the field mesh and the data is stored in the field data. \n\nParameters:\n\nT: type of data\nMT: type of mesh, e.g., Tuple{MeshType1, MeshType2, ...}\nN: number of dimensions\n\nMembers:\n\nmesh (MT): the mesh is a tuple of meshes.  The mesh should be an iterable object that contains an ordered list of grid points.   Examples are the \nMeshes defined in the MeshGrids module.\nUnitRange such as 1:10, etc.\nProduct of meshes MeshProduct defined in the MeshGrids module.\nIf a mesh is defined on a continuous manifold and supports the following methods, then one can perform interpolation, derivatives, etc. on the mesh:\nlocate(mesh, value): find the index of the closest grid point for given value;\nvolume(mesh, index): find the volume of grid space near the point at griven index.\nvolume(mesh, gridpoint): locate the corresponding index of a given grid point and than find the volume spanned by the grid point. \ndata (Array{T,N}): the data.\ndims: dimension of the data\n\n\n\n\n\n","category":"type"},{"location":"lib/mesharrays/#GreenFunc.MeshArrays.MeshArray-Tuple","page":"MeshArrays","title":"GreenFunc.MeshArrays.MeshArray","text":"function MeshArray(;\n    mesh...;\n    dtype = Float64,\n    data::Union{Nothing,AbstractArray}=nothing) where {T}\n\nCreate a Green struct. Its memeber dims is setted as the tuple consisting of the length of all meshes.\n\nArguments\n\nmesh: meshes of Green's function. See the docs of MeshArray for more details.  Mesh could be any iterable object, examples are vector, tuple, array, number, UnitRange (say, 1:5).\ndtype: data type of Green's function's value.\ndata: the data of the Green's function. By default, data is constructed to an unintialized Array with the dims size containing elements of dtype.\n\n\n\n\n\n","category":"method"},{"location":"lib/mesharrays/#GreenFunc.MeshArrays.MeshMatrix","page":"MeshArrays","title":"GreenFunc.MeshArrays.MeshMatrix","text":"MeshMatrix{T}\n\nAlias for MeshArray{T,2,MT}.\n\n\n\n\n\n","category":"type"},{"location":"lib/mesharrays/#GreenFunc.MeshArrays.MeshVector","page":"MeshArrays","title":"GreenFunc.MeshArrays.MeshVector","text":"MeshVector{T}\n\nAlias for MeshArray{T,1,MT}.\n\n\n\n\n\n","category":"type"},{"location":"lib/mesharrays/#Base.eltype-Union{Tuple{Type{GreenFunc.MeshArrays.AbstractMeshArray{T, N}}}, Tuple{N}, Tuple{T}} where {T, N}","page":"MeshArrays","title":"Base.eltype","text":"eltype(obj::AbstractMeshArray)\n\nReturn the type of the elements contained in obj.data.\n\n\n\n\n\n","category":"method"},{"location":"lib/mesharrays/#Base.getindex-Union{Tuple{N}, Tuple{MT}, Tuple{T}, Tuple{MeshArray{T, N, MT}, Vararg{Int64, N}}} where {T, MT, N}","page":"MeshArrays","title":"Base.getindex","text":"getindex(obj::MeshArray, inds...)\n\nReturn a subset of obj's data as specified by inds, where each inds may be, for example, an Int, an AbstractRange, or a Vector. \n\n\n\n\n\n","category":"method"},{"location":"lib/mesharrays/#Base.setindex!-Union{Tuple{N}, Tuple{MT}, Tuple{T}, Tuple{MeshArray{T, N, MT}, Any, Vararg{Int64, N}}} where {T, MT, N}","page":"MeshArrays","title":"Base.setindex!","text":"setindex!(obj::MeshArray, v, inds...)\nobj[inds...] = v\n\nStore values from array v within some subset of obj.data as specified by inds.\n\n\n\n\n\n","category":"method"},{"location":"lib/mesharrays/#Base.similar-Union{Tuple{S}, Tuple{N}, Tuple{MT}, Tuple{T}, Tuple{MeshArray{T, N, MT}, Type{S}}} where {T, MT, N, S}","page":"MeshArrays","title":"Base.similar","text":"Base.similar(obj::MeshArray{T,N,MT}, ::Type{S}) where {T,MT,N,S}\nBase.similar(obj::MeshArray{T,N,MT}) where {T,MT,N} = Base.similar(obj, T)\n\nReturn type:\n\nBase.similar(obj::MeshArray): Return a new MeshArray with the same meshes, and the uninitialized data of the same type as obj.data.\nBase.similar(obj::MeshArray, ::Type{S}): Return a new MeshArray with the same meshes, but with the uninitialized data of type S.\n\n\n\n\n\n","category":"method"},{"location":"lib/mesharrays/#Base.size-Tuple{GreenFunc.MeshArrays.AbstractMeshArray}","page":"MeshArrays","title":"Base.size","text":"size(obj::AbstractMeshArray)\n\nReturn a tuple containing the dimensions of obj.data (obj.dims).\n\n\n\n\n\n","category":"method"},{"location":"lib/mesharrays/#GreenFunc.MeshArrays._check-Tuple{MeshArray, MeshArray}","page":"MeshArrays","title":"GreenFunc.MeshArrays._check","text":"function _check(objL::MeshArray, objR::MeshArray)\n\nCheck if the Green's functions objL and objR are on the same meshes. Throw an AssertionError if any check is false.\n\n\n\n\n\n","category":"method"},{"location":"lib/mesharrays/#GreenFunc.MeshArrays.rank-Union{Tuple{Type{MeshArray{T, N, MT}}}, Tuple{N}, Tuple{MT}, Tuple{T}} where {T, MT, N}","page":"MeshArrays","title":"GreenFunc.MeshArrays.rank","text":"function rank(obj::MeshArray)\n\nReturn the dimension of obj.data (N).\n\n\n\n\n\n","category":"method"},{"location":"lib/greenfunc/#Green's-Functions","page":"GreenFunc","title":"Green's Functions","text":"","category":"section"},{"location":"lib/greenfunc/","page":"GreenFunc","title":"GreenFunc","text":"Modules = [GreenFunc]","category":"page"},{"location":"lib/greenfunc/#Base.:<<-Union{Tuple{MT2}, Tuple{MT1}, Tuple{N}, Tuple{T}, Tuple{MeshArray{T, N, MT1}, MeshArray{T, N, MT2}}} where {T, N, MT1, MT2}","page":"GreenFunc","title":"Base.:<<","text":"Base.:<<(objL::MeshArray, objR::MeshArray)\n\nDLR Fourier transform of functions on the first temporal grid (ImTime, ImFreq, or DLRFreq).\n\nIf objL and objR have identical temporal grid, objL<<objR assigns objR to objL.\nIf objL and objR have different temporal grid, one of them has to be in DLR space.\nIf objL is in DLR space, objL<<objR calculates the DLR spectral density of data in objR\nif objR is in DLR space, objL<<objR calculates the Green's function from the DLR spectral density in objR.\n\n\n\n\n\n","category":"method"},{"location":"lib/greenfunc/#GreenFunc.dlr_to_imtime-Union{Tuple{MeshArray{T, N, MT}}, Tuple{MT}, Tuple{N}, Tuple{T}, Tuple{MeshArray{T, N, MT}, Union{Nothing, AbstractVector}}} where {T, N, MT}","page":"GreenFunc","title":"GreenFunc.dlr_to_imtime","text":"function dlr_to_imfreq(mesharray[, tgrid; dim])\n\nTransform a Green's function in DLR to the imaginary-time domain. \n\nArguments\n\nmesharray: MeshArray in DLR space\ntgrid: The imaginary-time grid which the function transforms into. Default value is the imaginary-time grid from the DLRGrid from mesharray.mesh[dim].\ndim: The dimension of the temporal mesh. Default value is the first ImTime mesh.\n\n\n\n\n\n","category":"method"},{"location":"lib/greenfunc/#GreenFunc.imfreq_to_dlr-Union{Tuple{MeshArray{T, N, MT}}, Tuple{MT}, Tuple{N}, Tuple{T}} where {T, N, MT}","page":"GreenFunc","title":"GreenFunc.imfreq_to_dlr","text":"function imfreq_to_dlr(mesharray[; dim])\n\nCalculate the DLR spectral density of a Matsubara-frequency Green's function.\n\nArguments\n\nmesharray: MeshArray in the Matsubara-frequency domain.\ndim: The dimension of the mesh to be transformed. Default value is the first dimension with mesh type ImFreq.\n\n\n\n\n\n","category":"method"},{"location":"lib/greenfunc/#GreenFunc.imtime_to_dlr-Union{Tuple{MeshArray{T, N, MT}}, Tuple{MT}, Tuple{N}, Tuple{T}} where {T, N, MT}","page":"GreenFunc","title":"GreenFunc.imtime_to_dlr","text":"function imtime_to_dlr(mesharray[; dim])\n\nCalculate the DLR spectral density of an imaginary-time Green's function.\n\nArguments\n\nmesharray: MeshArray in the imaginary-time domain.\ndim: The dimension of the mesh to be transformed. Default value is the first dimension with mesh type ImTime.\n\n\n\n\n\n","category":"method"},{"location":"lib/greenfunc/#GreenFunc.to_dlr-Union{Tuple{MeshArray{T, N, MT}}, Tuple{MT}, Tuple{N}, Tuple{T}} where {T, N, MT}","page":"GreenFunc","title":"GreenFunc.to_dlr","text":"function to_dlr(mesharray[; dim])\n\nCalculate the DLR spectral density of an imaginary-time or Matsubara-frequency Green's function.\n\nArguments\n\nmesharray: MeshArray in the imaginary-time or the Matsubara-frequency domain.\ndim: The dimension of the mesh to be transformed. Default value is the first dimension with mesh type ImTime or ImFreq.\n\n\n\n\n\n","category":"method"},{"location":"lib/greenfunc/#GreenFunc.to_imfreq-Union{Tuple{MeshArray{T, N, MT}}, Tuple{MT}, Tuple{N}, Tuple{T}} where {T, N, MT}","page":"GreenFunc","title":"GreenFunc.to_imfreq","text":"function to_imfreq(mesharray[; dim])\n\nTransform a Green's function to the Matsubara-frequency domain.\n\nArguments\n\nmesharray: MeshArray in the imaginary-time, the Matsubara-frequency or the DLR frequency domain.\ndim: The dimension of the mesh to be transformed. Default value is the first dimension with mesh type DLRFreq, ImTime or ImFreq.\n\n\n\n\n\n","category":"method"},{"location":"lib/greenfunc/#GreenFunc.to_imtime-Union{Tuple{MeshArray{T, N, MT}}, Tuple{MT}, Tuple{N}, Tuple{T}} where {T, N, MT}","page":"GreenFunc","title":"GreenFunc.to_imtime","text":"function to_imtime(mesharray[; dim])\n\nTransform a Green's function to the imaginary-time domain.\n\nArguments\n\nmesharray: MeshArray in the imaginary-time, the Matsubara-frequency or the DLR frequency domain.\ndim: The dimension of the mesh to be transformed. Default value is the first dimension with mesh type DLRFreq, ImTime or ImFreq.\n\n\n\n\n\n","category":"method"},{"location":"lib/meshgrids/#MeshGrids","page":"MeshGrids","title":"MeshGrids","text":"","category":"section"},{"location":"lib/meshgrids/","page":"MeshGrids","title":"MeshGrids","text":"Modules = [GreenFunc.MeshGrids]","category":"page"},{"location":"lib/meshgrids/#GreenFunc.MeshGrids.DLRFreq","page":"MeshGrids","title":"GreenFunc.MeshGrids.DLRFreq","text":"function DLRFreq(β, isFermi::Bool=false;\n    dtype=Float64,\n    rtol=1e-12,\n    Euv=1000 / β,\n    sym=:none,\n    rebuild=false,\n    dlr::Union{DLRGrid,Nothing}=nothing\n)\n\nCreate a DLRFreq struct from parameters.\n\nArguments\n\nβ: inverse temperature.\nisFermi: the statistics for particles is fermionic or not. False by default.\ndtype: type of β and Euv.\nrtol: tolerance absolute error. By default, rtol = 1e-12.\nEuv: the UV energy scale of the spectral density. By default, Euv = 1000 / β.\nsymmetry: :ph for particle-hole symmetric, :pha for particle-hole symmetry, and :none for no symmetry. By default, sym = :none.\nrebuild: if no dlr is input, set false to load DLRGrid from the file; set true to recalculate the DLRGrid on the fly. By default, rebuild = false.\n\n\n\n\n\n","category":"type"},{"location":"lib/meshgrids/#GreenFunc.MeshGrids.DLRFreq-2","page":"MeshGrids","title":"GreenFunc.MeshGrids.DLRFreq","text":"struct DLRFreq{T<:Real} <: TemporalGrid{Int}\n\nDiscrete-Lehmann-representation grid for Green's functions. \n\nParameters\n\nT: type of the grid point, β and Euv.\n\nMembers\n\ndlr: built-in DLR grid.\ngrid: 1D grid of time axis, with locate, volume, and AbstractArray interface implemented. It should be grid of Int for ImFreq, and DLRGrid for DLRFreq.\nβ: inverse temperature.\nEuv:  the UV energy scale of the spectral density.\nrtol: tolerance absolute error.\nsymmetry: :ph for particle-hole symmetric, :pha for particle-hole symmetry, and :none for no symmetry. By default, sym = :none.\nisFermi: the statistics for particles. \n\n\n\n\n\n","category":"type"},{"location":"lib/meshgrids/#GreenFunc.MeshGrids.DLRFreq-Tuple{Lehmann.DLRGrid}","page":"MeshGrids","title":"GreenFunc.MeshGrids.DLRFreq","text":"function DLRFreq(dlr::DLRGrid)\n\nCreate a DLRFreq struct from DLRGrid.\n\nArguments\n\ndlr: 1D DLR grid.\n\n\n\n\n\n","category":"method"},{"location":"lib/meshgrids/#GreenFunc.MeshGrids.ImFreq","page":"MeshGrids","title":"GreenFunc.MeshGrids.ImFreq","text":"function ImFreq(β, isFermi::Bool=false;\n    dtype=Float64,\n    Euv=1000 / β,\n    rtol=1e-12,\n    symmetry=:none,\n    grid::Union{AbstractGrid,AbstractVector,Nothing}=nothing\n)\n\nCreate a ImFreq struct from parameters.\n\nArguments\n\nβ: inverse temperature.\nisFermi: the statistics for particles is fermionic or not. False by default.\ndtype: type of β and Euv. By default, dtype = Float64.\nEuv: the UV energy scale of the spectral density. By default, Euv = 1000 / β.\nsymmetry: :ph for particle-hole symmetric, :pha for particle-hole symmetry, and :none for no symmetry. By default, sym = :none.\ngrid: 1D Matsubara-frequency integer-valued grid as a AbstractVector or CompositeGrids.AbstractGrid. By default, a optimized grid built in DLR is used.\n\n\n\n\n\n","category":"type"},{"location":"lib/meshgrids/#GreenFunc.MeshGrids.ImFreq-2","page":"MeshGrids","title":"GreenFunc.MeshGrids.ImFreq","text":"struct ImFreq{T, G, R} <: TemporalGrid{Int}\n\nImaginary-frequency grid for Green's functions. \n\nParameters\n\nT<:Real: type of the grid point, β and Euv.\nG<:AbstractGrid{T}: type of 1D grid with T as the grid point type.\nR: type of the representation.\nREV: access the grid in reverse order or not.\n\nMembers\n\ngrid: 1D grid of time axis, with locate, volume, and AbstractArray interface implemented. Always in ascend order It should be grid of Int for ImFreq, and DLRGrid for DLRFreq.\nβ: inverse temperature.\nEuv:  the UV energy scale of the spectral density.\nisFermi: the statistics for particles is fermionic or not.\nsymmetry: :ph for particle-hole symmetric, :pha for particle-hole symmetry, and :none for no symmetry. By default, sym = :none.\nrtol: relative tolerance\nrepresentation: the representation of the Green's function.\n\n\n\n\n\n","category":"type"},{"location":"lib/meshgrids/#GreenFunc.MeshGrids.ImFreq-Tuple{Lehmann.DLRGrid}","page":"MeshGrids","title":"GreenFunc.MeshGrids.ImFreq","text":"function ImFreq(dlr::DLRGrid;\n    dtype=Float64,\n    grid::Union{AbstractGrid,AbstractVector}=SimpleG.Arbitrary{Int}(dlr.n)\n)\n\nConstruct ImFreq from a DLRGrid, with a given grid. By default, grid is the Matsubara-frequency points from DLRGrid.\n\n\n\n\n\n","category":"method"},{"location":"lib/meshgrids/#GreenFunc.MeshGrids.ImTime","page":"MeshGrids","title":"GreenFunc.MeshGrids.ImTime","text":"function ImTime(β, isFermi::Bool=false;\n    dtype=Float64,\n    rtol=1e-12,\n    Euv=1000 / β,\n    symmetry=:none,\n    grid::Union{AbstractGrid,AbstractVector,Nothing}=nothing\n)\n\nCreate a ImTime struct.\n\nArguments\n\nβ: inverse temperature.\nisFermi: the statistics for particles is fermionic or not. False by default.\ndtype: type of the grid point. By default, dtype = Float64.\nEuv: the UV energy scale of the spectral density. By default, Euv = 1000 / β.\nsymmetry: :ph for particle-hole symmetric, :pha for particle-hole symmetry, and :none for no symmetry. By default, sym = :none.\ngrid: 1D time grid as a AbstractVector or CompositeGrids.AbstractGrid. By default, a optimized grid built in DLR is used.\n\n\n\n\n\n","category":"type"},{"location":"lib/meshgrids/#GreenFunc.MeshGrids.ImTime-2","page":"MeshGrids","title":"GreenFunc.MeshGrids.ImTime","text":"struct ImTime{T, G, R} <: TemporalGrid{T}\n\nTime grid for Green's functions.\n\nParameters\n\nT<:Real: type of the grid point, β and Euv.\nG<:AbstractGrid{T}: type of 1D grid with T as the grid point type.\nR: type of the representation.\nREV: access the grid in reverse order or not.\n\nMembers\n\ngrid: 1D grid of time axis, with locate, volume, and AbstractArray interface implemented. Always in ascend order. It should be grid of Int for ImFreq, and DLRGrid for DLRFreq.\nβ: inverse temperature.\nEuv:  the UV energy scale of the spectral density.\nisFermi: the statistics for particles is fermionic or not.\nsymmetry: :ph for particle-hole symmetric, :pha for particle-hole symmetry, and :none for no symmetry. By default, sym = :none.\nrtol: relative tolerance\nrepresentation: the representation of the Green's function.\n\n\n\n\n\n","category":"type"},{"location":"lib/meshgrids/#GreenFunc.MeshGrids.ImTime-Tuple{Lehmann.DLRGrid}","page":"MeshGrids","title":"GreenFunc.MeshGrids.ImTime","text":"function ImTime(dlr::DLRGrid;\n                dtype=Float64,\n                grid::Union{AbstractGrid,AbstractVector}=SimpleG.Arbitrary{dtype}(dlr.τ))\n\nConstruct ImTime from a DLRGrid, with a given grid. By default, grid is the imaginary-time grid points from DLRGrid.\n\n\n\n\n\n","category":"method"},{"location":"lib/meshgrids/#GreenFunc.MeshGrids.MeshProduct","page":"MeshGrids","title":"GreenFunc.MeshGrids.MeshProduct","text":"MeshProduct{MT,N}\nMeshProduct(vargs...)\n\nThe Cartesian Mesh product:\n\nParameters:\n\nMT: Type of meshes\nN : Number of meshes\n\nMembers:\n\nmeshes: The list of Meshes in the MeshProduct\ndims: A tuple of the length of the mesh factors\n\n\n\n\n\n","category":"type"},{"location":"lib/meshgrids/#Base.floor-Union{Tuple{T}, Tuple{TemporalGrid{T, false}, Any}} where T","page":"MeshGrids","title":"Base.floor","text":"Base.floor(tg::TemporalGrid{T,false}, pos) where {T} = floor(tg.grid, pos) #TODO: how to implement?\nBase.floor(tg::TemporalGrid{T,true}, pos) where {T} = length(tg) - floor(tg.grid, pos) #TODO: how to implement?\n\nIf the grid is in ascend order, then floor returns the largest index that the grid point is smaller than pos. If the grid is in descend order, then floor returns the largest index that the grid point is larger than pos.\n\nIn both cases, the returned index is in the range [1, length(tg)-1]\n\n\n\n\n\n","category":"method"},{"location":"lib/meshgrids/#Base.getindex-Union{Tuple{R}, Tuple{G}, Tuple{T}, Tuple{ImFreq{T, G, R, false}, Int64}} where {T, G, R}","page":"MeshGrids","title":"Base.getindex","text":"getindex(g::ImFreq{T, G, R, REV}, I::Int)\n\nEquivalent to g[I], get the real-valued Matsubara frequency of the Ith point in the grid. For fermion, return (2g[I]+1)π/β, for boson, return 2g[I]*π/β.\n\nIf REV = true, then index in the reversed order, namely I will be replaced with length(g) - I + 1.\n\nIf you need the integer-valued frequency, use g.grid[I] instead.\n\n\n\n\n\n","category":"method"},{"location":"lib/meshgrids/#Base.length-Tuple{MeshProduct}","page":"MeshGrids","title":"Base.length","text":"function Base.length(obj::MeshProduct)\n\nReturn the number of grids of the MeshProduct.\n\n\n\n\n\n","category":"method"},{"location":"lib/meshgrids/#Base.show-Tuple{IO, DLRFreq}","page":"MeshGrids","title":"Base.show","text":"show(io::IO, tg::DLRFreq)\n\nWrite a text representation of the DLR grid tg to the output stream io.\n\n\n\n\n\n","category":"method"},{"location":"lib/meshgrids/#Base.show-Tuple{IO, ImFreq}","page":"MeshGrids","title":"Base.show","text":"show(io::IO, tg::ImFreq)\n\nWrite a text representation of the Imaginary-frequency grid tg to the output stream io.\n\n\n\n\n\n","category":"method"},{"location":"lib/meshgrids/#Base.show-Tuple{IO, ImTime}","page":"MeshGrids","title":"Base.show","text":"show(io::IO, tg::ImTime)\n\nWrite a text representation of the Imaginary-time grid tg to the output stream io.\n\n\n\n\n\n","category":"method"},{"location":"lib/meshgrids/#Base.show-Tuple{IO, MeshProduct}","page":"MeshGrids","title":"Base.show","text":"function Base.show(io::IO, obj::MeshProduct)\n\nPrint the MeshProduct.\n\n\n\n\n\n","category":"method"},{"location":"lib/meshgrids/#Base.size-Tuple{MeshProduct, Int64}","page":"MeshGrids","title":"Base.size","text":"function Base.size(obj::MeshProduct, I::Int)\n\nReturn the length of the specific Ith mesh factor of the MeshProduct.\n\n\n\n\n\n","category":"method"},{"location":"lib/meshgrids/#Base.size-Tuple{MeshProduct}","page":"MeshGrids","title":"Base.size","text":"function Base.size(obj::MeshProduct, I::Int)\n\nReturn the length of the specific Ith mesh factor of the MeshProduct.\n\n\n\n\n\n","category":"method"},{"location":"lib/meshgrids/#GreenFunc.MeshGrids.index_to_linear-Union{Tuple{N}, Tuple{MT}, Tuple{MeshProduct{MT, N}, Vararg{Any}}} where {MT, N}","page":"MeshGrids","title":"GreenFunc.MeshGrids.index_to_linear","text":"function index_to_linear(obj::MeshProduct, index...)\n\nConvert a tuple of the indexes of each mesh to a single linear index of the MeshProduct.\n\nArgument:\n\nobj: The MeshProduct object\nindex...: N indexes of the mesh factor, where N is the number of mesh factor\n\n\n\n\n\n","category":"method"},{"location":"lib/meshgrids/#GreenFunc.MeshGrids.linear_to_index-Union{Tuple{N}, Tuple{MT}, Tuple{MeshProduct{MT, N}, Int64}} where {MT, N}","page":"MeshGrids","title":"GreenFunc.MeshGrids.linear_to_index","text":"function linear_to_index(obj::MeshProduct, I::Int)\n\nConvert the single linear index of the MeshProduct to a tuple of indexes of each mesh. \n\nArgument:\n\nobj: The MeshProduct object\nI: The linear index of the MeshProduct\n\n\n\n\n\n","category":"method"},{"location":"lib/meshgrids/#GreenFunc.MeshGrids.locate-Tuple{ImFreq, Int64}","page":"MeshGrids","title":"GreenFunc.MeshGrids.locate","text":"locate(tg::ImFreq, n::Int)\nlocate(tg::ImFreq, ωn)\n\nFind the location in tg.grid for the Matsubara frequency ωn or the integer n.\n\n\n\n\n\n","category":"method"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = GreenFunc","category":"page"},{"location":"#GreenFunc","page":"Home","title":"GreenFunc","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"(Image: Stable) (Image: Dev) (Image: Build Status) (Image: Coverage)","category":"page"},{"location":"","page":"Home","title":"Home","text":"GreenFunc.jl is a differentiable numerical framework to manipulate multidimensional Green's functions.","category":"page"},{"location":"#Features","page":"Home","title":"Features","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"MeshArray type as an array defined on meshes, which provides a generic data structure for Green's functions, vertex functions or any other correlation/response functions.\nStructured (non-)uniform Brillouin Zone meshes powered by the package BrillouinZoneMeshes.jl.\nStructured (non-)uniform temporal meshes for (imaginary-)time or (Matsubara-)frequency domains powered by the pacakge CompositeGrids.jl.\nCompat representation based on the Discrete Lehmann representation (DLR) powered by the package Lehmann.jl.\nAccurate and fast Fourier transform.\nInterface to the TRIQS library.","category":"page"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"This package has been registered. So, simply type","category":"page"},{"location":"","page":"Home","title":"Home","text":"import Pkg; Pkg.add(\"GreenFunc\")","category":"page"},{"location":"","page":"Home","title":"Home","text":"in the Julia REPL to install.","category":"page"},{"location":"#Basic-Usage","page":"Home","title":"Basic Usage","text":"","category":"section"},{"location":"#Example-1:-Green's-function-of-a-single-level","page":"Home","title":"Example 1: Green's function of a single level","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"We first show how to use MeshArray to present Green's function of a single-level quantum system filled with spinless fermionic particles. We assume that the system could exchange particles and energy with the environment so that it's equilibrium state is a grand canonical ensemble. The single-particle Green's function then has a simple form in Matsubara-frequency representation:  G(ωₙ) = frac1(iωₙ - E) where E is the level energy. We show how to generate and manipulate this Green's function.","category":"page"},{"location":"","page":"Home","title":"Home","text":"using GreenFunc\n\n# inverse temperature and the level energy\nβ = 100.0; E = 1.0\n# UV energy cutoff is 100 times larger than the level energy\nωₙ_mesh = MeshGrids.ImFreq(100.0, FERMION; Euv = 100E)\n# Green's function defined on the ωₙ_mesh\nGn = MeshArray(ωₙ_mesh; dtype=ComplexF64)\n\nfor (n, ωₙ) in enumerate(Gn.mesh[1])\n    Gn[n] = 1/(ωₙ*im - E)\nend","category":"page"},{"location":"","page":"Home","title":"Home","text":"Green's function describes correlations between two or more spacetime events. The spacetime continuum needs to be discretized into spatial and temporal meshes. This example demonstrates how to define a one-body Green's function on a temporal mesh. The package provides three types of temporal meshes: imaginary-time grid, Matsubara-frequency grid, and DLR grid. The latter provides a generic compressed representation for Green's functions (We will show how to use DLR later).  Correspondingly, They can be created with the ImTime, ImFreq, and DLRFreq methods. The user needs to specify the inverse temperature, whether the particle is fermion or boson (using the constant FERMION or BOSON). Internally, a set of non-uniform grid points optimized for the given inverse temperature and the cutoff energy will be created with the given parameters.\nOnce the meshes are created, one can define a MeshArray on them to represent the Green's function Gn. The constructor of MeshArray takes a set of meshes and initializes a multi-dimensional array. Each mesh corresponds to one dimension of the array. The data type of the MeshArray is specified by the optional keyword argument dtype, which is set to Float64 by default. You can access the meshes (stored as a tuple) with Gn.mesh and the array data with Gn.data.\nBy default, Gn.data is left undefined if not specified by the user. To initialize it, one can either use the optional keyword argument data in the constructor or use the iterator interface of the meshes and the MeshArray. ","category":"page"},{"location":"#Example-2:-Green's-function-of-a-free-electron-gas","page":"Home","title":"Example 2: Green's function of a free electron gas","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Now let us show how to create a Green's function of a free electron gas. Unlike the spinless fermionic particle, the electron is a spin-1/2 particle so that it has two inner states. In free space, it has a kinetic energy ϵ_q = q^2-E (we use the unit where m_e = 12). The Green's function in Matsubara-frequency space is then given by the following equation: G_n = G_sigma_1 sigma_2(qomega_n) = frac1i omega_n - epsilon_q, where sigma_i denotes the spins of the incoming and the outgoing electron in the propagator. We inherit the Matsubara-frequency grid from the first example. We show how to use the CompositeGrids package to generate momentum grids and how to treat the multiple inner states and the meshes with MeshArray.","category":"page"},{"location":"","page":"Home","title":"Home","text":"using GreenFunc, CompositeGrids\n\n# inverse temperature and the level energy\nβ = 100.0; E = 1.0\n# UV energy cutoff is 100 times larger than the level energy\nωₙ_mesh = MeshGrids.ImFreq(100.0, FERMION; Euv = 100E)\n# initialze an uniform momentum grid\nkmesh = SimpleGrid.Uniform{Float64}([0.0, 10.0], 50)\n# Green's function of free electron gas with 2x2 innerstates\nG_n =  MeshArray(1:2, 1:2, kmesh, ωₙ_mesh; dtype=ComplexF64)\n\nfor ind in eachindex(G_n)\n    q = G_n.mesh[3][ind[3]]\n    ω_n = G_n.mesh[4][ind[4]]\n    G_n[ind] = 1/(ω_n*im - (q^2-E))\nend","category":"page"},{"location":"","page":"Home","title":"Home","text":"One can generate various types of grids with the CompositeGrids package. The SimpleGrid module here provides several basic grids, such as uniform grids and logarithmically dense grids. TheUniform method here generates a 1D linearly spaced grid. The user has to specify the number of grid points N and the boundary points [min, max]. One can also combine arbitrary numbers of SimpleGrid subgrids with a user-specified pattern defined by a panel grid. These more advanced grids optimized for different purposes can be found in this link.\nThe constructor of MeshArray can take any iterable objects as one of its meshes. Therefore for discrete inner states such as spins, one can simply use a 1:2, which is a UnitRange{Int64} object.","category":"page"},{"location":"#Example-3:-Green's-function-of-a-Hubbard-lattice","page":"Home","title":"Example 3: Green's function of a Hubbard lattice","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Now we show how to generate a multi-dimensional Green's function on a Brillouin Zone meshe. We calculate the Green's function of a free spinless Fermi gas on a square lattice. It has a tight-binding dispersion epsilon_q = -2t(cos(q_x)+cos(q_y)), which gives G(q omega_n) = frac1iomega_n - epsilon_q. The momentum is defined on the first Brillouin zone captured by a 2D k-mesh.","category":"page"},{"location":"","page":"Home","title":"Home","text":"using GreenFunc\nusing GreenFunc: BrillouinZoneMeshes\n\nDIM, nk = 2, 8\nlatvec = [1.0 0.0; 0.0 1.0] .* 2π\nbzmesh = BrillouinZoneMeshes.BaseMesh.UniformMesh{DIM, nk}([0.0, 0.0], latvec)\nωₙmesh = ImFreq(10.0, FERMION)\ng_freq =  MeshArray(bzmesh, ωₙmesh; dtype=ComplexF64)\n\nt = 1.0\nfor ind in eachindex(g_freq)\n    q = g_freq.mesh[1][ind[1]]\n    ωₙ = g_freq.mesh[2][ind[2]]\n    g_freq[ind] = 1/(ωₙ*im - (-2*t*sum(cos.(q))))\nend","category":"page"},{"location":"","page":"Home","title":"Home","text":"For lattice systems with multi-dimensional Brillouin zone, the momentum grids internally generated with the BrillouinZoneMeshes.jl package. Here a UniformMesh{DIM,N}(origin, latvec) generates a linearly spaced momentum mesh on the first Brillouin zone defined by origin and lattice vectors given. For more detail, see the link.","category":"page"},{"location":"#Example-4:-Fourier-Transform-of-Green's-function-with-DLR","page":"Home","title":"Example 4:  Fourier Transform of Green's function with DLR","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"DLR provides a compact representation for one-body Green's functions. At a temperature T and an accuracy level epsilon, it represents a generic Green's function with only log (1T) log (1epsilon) basis functions labeled by a set of real frequency grid points. It enables fast Fourier transform and interpolation between the imaginary-time and the Matsubara-frequency representations with a cost O(log (1T) log (1epsilon)). GreenFunc.jl provide DLR through the package Lehmann.jl.","category":"page"},{"location":"","page":"Home","title":"Home","text":"In the following example, we demonstrate how to perform DLR-based Fourier transform in GreenFunc.jl between the imaginary-time and the Matsubara-frequency domains back and forth through the DLR representation.","category":"page"},{"location":"","page":"Home","title":"Home","text":"using GreenFunc, CompositeGrids\n\nβ = 100.0; E = 1.0 # inverse temperature and the level energy\nωₙ_mesh = ImFreq(100.0, FERMION; Euv = 100E) # UV energy cutoff is 100 times larger than the level energy\nkmesh = SimpleGrid.Uniform{Float64}([0.0, 10.0], 50); # initialze an uniform momentum grid\nG_n =  MeshArray(1:2, 1:2, kmesh, ωₙ_mesh; dtype=ComplexF64); # Green's function of free electron gas with 2x2 innerstates\n\nfor ind in eachindex(G_n)\n    q = G_n.mesh[3][ind[3]]\n    ω_n = G_n.mesh[4][ind[4]]\n    G_n[ind] = 1/(im*ω_n - (q^2-E))\nend\n\nG_dlr = to_dlr(G_n) # convert G_n to DLR space\nG_tau = to_imtime(G_dlr) # convert G_dlr to the imaginary-time domain\n\n#alternative, you can use the pipe operator\nG_tau = G_n |> to_dlr |> to_imtime #Fourier transform to (k, tau) domain\n","category":"page"},{"location":"","page":"Home","title":"Home","text":"The imaginary-time Green's function after the Fourier transform shoud be consistent with the analytic solution G_tau = -e^-tau epsilon_q(1+e^-beta epsilon_q).","category":"page"},{"location":"","page":"Home","title":"Home","text":"For any Green's function that has at least one imaginary-temporal grid (ImTime, ImFreq, and DLRFreq) in meshes, we provide a set of operations (to_dlr, to_imfreq and to_imtime) to bridge the DLR space with imaginary-time and Matsubara-frequency space. By default, all these functions find the dimension of the imaginary-temporal mesh within Green's function meshes and perform the transformation with respect to it. Alternatively, one can specify the dimension with the optional keyword argument dim. Be careful that the original version of DLR is only guaranteed to work with one-body Green's function.\nOnce a spectral density G_dlr in DLR space is obtained, one can use to_imfreq or to_imtime methods to reconstruct the Green's function in the corresponding space. By default, to_imfreq and to_imtime uses an optimized imaginary-time or Matsubara-frequency grid from the DLR. User can assign a target imaginary-time or Matsubara-frequency grid if necessary.   \nCombining to_dlr, to_imfreq and to_imtime allows both interpolation as well as Fourier transform.\nSince the spectral density G_dlr can be reused whenever the user wants to change the grid points of Green's function (normally through interpolation that lost more accuracy than DLR transform), we encourage the user always to keep the G_dlr objects. If the intermediate DLR Green's function is not needed, the user can use piping operator |> as shown to do Fourier transform directly between ImFreq and ImTime in one line.","category":"page"},{"location":"#Interface-with-TRIQS","page":"Home","title":"Interface with TRIQS","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"TRIQS (Toolbox for Research on Interacting Quantum Systems) is a scientific project providing a set of C++ and Python libraries for the study of interacting quantum systems. We provide a direct interface to convert TRIQS objects, such as the temporal meshes, the Brillouin zone meshes, and the  multi-dimensional (blocked) Green's functions, to the equivalent objects in our package. It would help TRIQS users to make use of our package without worrying about the different internal data structures.","category":"page"},{"location":"","page":"Home","title":"Home","text":"The interface is provided by an independent package NEFTInterface.jl. We provide several examples of interfacing TRIQS and GreenFunc.jl in the NEFTInterface.jl README.","category":"page"}]
}
