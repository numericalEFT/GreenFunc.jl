var documenterSearchIndex = {"docs":
[{"location":"lib/deprecated/#Deprecated-Green's-Functions-API","page":"Deprecated","title":"Deprecated Green's Functions API","text":"","category":"section"},{"location":"lib/deprecated/","page":"Deprecated","title":"Deprecated","text":"Modules = [GreenFunc.Deprecated]","category":"page"},{"location":"lib/deprecated/#GreenFunc.Deprecated.Green2DLR","page":"Deprecated","title":"GreenFunc.Deprecated.Green2DLR","text":"Green's function with two external legs that has in-built Discrete Lehmann Representation. #Parameters:\n\n'T': type of data\n'TType': type of time domain, TType<:TimeDomain\n'TGT': type of time grid\n'SGT': type of space grid\n\n#Members:\n\n'name': Name of green's function\n'color': Number of different species of Green's function (such as different spin values)\n'dlrGrid': In-built Discrete Lehmann Representation\n'timeGrid': Time or Frequency grid\n'spaceType': Whether the Green's function is in coordinate space/momentum space\n'spaceGrid': Coordinate or momentum grid\n'instant': Instantaneous part of Green's function that is proportional to δ(τ) in τ space.\n'dynamic': Dynamic part of Green's function\n'instantError': Error of instantaneous part\n'dynamicError': Error of dynamic part\n\n\n\n\n\n","category":"type"},{"location":"lib/deprecated/#GreenFunc.Deprecated.GreenSym2DLR","page":"Deprecated","title":"GreenFunc.Deprecated.GreenSym2DLR","text":"Symmetrized Green's function with two external legs that has in-built Discrete Lehmann Representation. The real and imaginary parts are saved separately on corresponding symmetrized DLR grids.  #Parameters:\n\n'T': type of data\n'TType': type of time domain, TType<:TimeDomain\n'TGT': type of time grid\n'SGT': type of space grid\n\n#Members:\n\n'name': Name of green's function\n'color': Number of different species of Green's function (such as different spin values)\n'dlrGrid': In-built Discrete Lehmann Representation\n'timeGrid': Time or Frequency grid\n'spaceType': Whether the Green's function is in coordinate space/momentum space\n'spaceGrid': Coordinate or momentum grid\n'instant': Instantaneous part of Green's function that is proportional to δ(τ) in τ space.\n'dynamic': Dynamic part of Green's function\n'instantError': Error of instantaneous part\n'dynamicError': Error of dynamic part\n\n\n\n\n\n","category":"type"},{"location":"lib/deprecated/#GreenFunc.Deprecated._dynamic-Union{Tuple{SGT}, Tuple{TGT}, Tuple{TT}, Tuple{DT}, Tuple{Any, Any}} where {DT, TT, TGT<:CompositeGrids.SimpleG.AbstractGrid, SGT<:CompositeGrids.SimpleG.AbstractGrid}","page":"Deprecated","title":"GreenFunc.Deprecated._dynamic","text":"function dynamic(green::Union{Green2DLR{DT,TT,TGT,SGT},GreenSym2DLR{DT,TT,TGT,SGT}}, time, space, color1::Int, color2::Int, timeMethod::TM , spaceMethod::SM) where {DT,TT,TGT<:CompositeGrids.AbstractGrid,SGT<:CompositeGrids.AbstractGrid,TM,SM}\n\nFind value of Green's function's dynamic part at given color and k/x by interpolation. Interpolation method is by default depending on the grid, but could also be chosen to be linear.\n\n#Argument\n\n'green': Green's function\n'time': Target τ/ω_n point\n'space': Target k/x point\n'color1': Target color1\n'color2': Target color2\n'timeMethod': Method of interpolation for time\n'spaceMethod': Method of interpolation for space \n\n\n\n\n\n","category":"method"},{"location":"lib/deprecated/#GreenFunc.Deprecated.instant-Union{Tuple{}, Tuple{SM}, Tuple{SGT}, Tuple{TGT}, Tuple{TT}, Tuple{DT}} where {DT, TT, TGT, SGT, SM}","page":"Deprecated","title":"GreenFunc.Deprecated.instant","text":"function instant(green::Green2DLR{DT,TT,TGT,SGT}, space, color1::Int, color2::Int=color1; spaceMethod::SM = DEFAULTINTERP) where {DT,TT,TGT,SGT,SM}\n\nFind value of Green's function's instant part at given color and k/x by interpolation. Interpolation method is by default depending on the grid, but could also be chosen to be linear.\n\n#Argument\n\n'green': Green's function\n'space': Target k/x point\n'color1': Target color1\n'color2': Target color2\n'spaceMethod': Method of interpolation for space. \n\n\n\n\n\n","category":"method"},{"location":"lib/deprecated/#GreenFunc.Deprecated.toDLR-Tuple{GreenFunc.Deprecated.Green2DLR}","page":"Deprecated","title":"GreenFunc.Deprecated.toDLR","text":"function toDLR(green::Green2DLR)\n\nConvert Green's function to dlr space.\n\n#Arguements\n\n'green': Original Green's function\n\n\n\n\n\n","category":"method"},{"location":"lib/deprecated/#GreenFunc.Deprecated.toDLR-Tuple{GreenFunc.Deprecated.GreenSym2DLR}","page":"Deprecated","title":"GreenFunc.Deprecated.toDLR","text":"function toDLR(green::Green2DLR)\n\nConvert Green's function to dlr space.\n\n#Arguements\n\n'green': Original Green's function\n\n\n\n\n\n","category":"method"},{"location":"lib/deprecated/#GreenFunc.Deprecated.toMatFreq","page":"Deprecated","title":"GreenFunc.Deprecated.toMatFreq","text":"function toMatFreq(green::Green2DLR, targetGrid =  green.dlrGrid.n)\n\nConvert Green's function to matfreq space by Fourier transform. If green is already in matfreq space then it will be interpolated to the new grid.\n\n#Arguements\n\n'green': Original Green's function\n'targetGrid': Grid of outcome Green's function. Default: DLR n grid\n\n\n\n\n\n","category":"function"},{"location":"lib/deprecated/#GreenFunc.Deprecated.toMatFreq-2","page":"Deprecated","title":"GreenFunc.Deprecated.toMatFreq","text":"function toMatFreq(green::Green2DLR, targetGrid =  green.dlrGrid.n)\n\nConvert Green's function to matfreq space by Fourier transform. If green is already in matfreq space then it will be interpolated to the new grid.\n\n#Arguements\n\n'green': Original Green's function\n'targetGrid': Grid of outcome Green's function. Default: DLR n grid\n\n\n\n\n\n","category":"function"},{"location":"lib/deprecated/#GreenFunc.Deprecated.toTau","page":"Deprecated","title":"GreenFunc.Deprecated.toTau","text":"function toTau(green::Green2DLR, targetGrid =  green.dlrGrid.τ)\n\nConvert Green's function to τ space by Fourier transform. If green is already in τ space then it will be interpolated to the new grid.\n\n#Arguements\n\n'green': Original Green's function\n'targetGrid': Grid of outcome Green's function. Default: DLR τ grid\n\n\n\n\n\n","category":"function"},{"location":"lib/deprecated/#GreenFunc.Deprecated.toTau-2","page":"Deprecated","title":"GreenFunc.Deprecated.toTau","text":"function toTau(green::Green2DLR, targetGrid =  green.dlrGrid.τ)\n\nConvert Green's function to τ space by Fourier transform. If green is already in τ space then it will be interpolated to the new grid.\n\n#Arguements\n\n'green': Original Green's function\n'targetGrid': Grid of outcome Green's function. Default: DLR τ grid\n\n\n\n\n\n","category":"function"},{"location":"lib/mesharrays/#Array-with-MeshGrids","page":"MeshArrays","title":"Array with MeshGrids","text":"","category":"section"},{"location":"lib/mesharrays/","page":"MeshArrays","title":"MeshArrays","text":"Modules = [GreenFunc.MeshArrays]","category":"page"},{"location":"lib/mesharrays/#GreenFunc.MeshArrays.MeshArray","page":"MeshArrays","title":"GreenFunc.MeshArrays.MeshArray","text":"struct MeshArray{T,N,MT} <: AbstractMeshArray{T,N}\n\nMulti-dimensional array that is defined on a mesh.  The mesh is a tuple of meshgrid objects.  The mesh is stored in the field mesh and the data is stored in the field data. \n\nParameters:\n\nT: type of data\nMT: type of mesh, e.g., Tuple{MeshType1, MeshType2, ...}\nN: number of dimensions\n\nMembers:\n\nmesh (MT): the mesh is a tuple of meshes.    The mesh should be an iterable object that contains an ordered list of grid points.   Examples are the \nMeshes defined in the MeshGrids module.\nUnitRange such as 1:10, etc.\nProduct of meshes MeshProduct defined in the MeshGrids module.\nIf a mesh is defined on a continuous manifold and supports the following methods, then one can perform interpolation, derivatives, etc. on the mesh:\nlocate(mesh, value): find the index of the closest grid point for given value;\nvolume(mesh, index): find the volume of grid space near the point at griven index.\nvolume(mesh, gridpoint): locate the corresponding index of a given grid point and than find the volume spanned by the grid point. \ndata: Array{T,N}: the data.\ndims: dimension of the data\n\n\n\n\n\n","category":"type"},{"location":"lib/mesharrays/#GreenFunc.MeshArrays.MeshArray-Tuple","page":"MeshArrays","title":"GreenFunc.MeshArrays.MeshArray","text":"function MeshArray(;\n    mesh...;\n    dtype = Float64,\n    data::Union{Nothing,AbstractArray}=nothing) where {T}\n\nCreate a Green struct. Its memeber dims is setted as the tuple consisting of the length of all meshes.\n\nArguments\n\nmesh: meshes of Green's function. See the docs of MeshArray for more details.  Mesh could be any iterable object, examples are vector, tuple, array, number, UnitRange (say, 1:5).\ndtype: data type of Green's function's value.\ndata: the data of the Green's function. By default, data is constructed to an unintialized Array with the dims size containing elements of dtype.\n\n\n\n\n\n","category":"method"},{"location":"lib/mesharrays/#GreenFunc.MeshArrays.MeshMatrix","page":"MeshArrays","title":"GreenFunc.MeshArrays.MeshMatrix","text":"MeshMatrix{T}\n\nAlias for MeshArray{T,2,MT}.\n\n\n\n\n\n","category":"type"},{"location":"lib/mesharrays/#GreenFunc.MeshArrays.MeshVector","page":"MeshArrays","title":"GreenFunc.MeshArrays.MeshVector","text":"MeshVector{T}\n\nAlias for MeshArray{T,1,MT}.\n\n\n\n\n\n","category":"type"},{"location":"lib/mesharrays/#Base.eltype-Union{Tuple{Type{GreenFunc.MeshArrays.AbstractMeshArray{T, N}}}, Tuple{N}, Tuple{T}} where {T, N}","page":"MeshArrays","title":"Base.eltype","text":"eltype(obj::AbstractMeshArray)\n\nReturn the type of the elements contained in obj.data.\n\n\n\n\n\n","category":"method"},{"location":"lib/mesharrays/#Base.getindex-Union{Tuple{N}, Tuple{MT}, Tuple{T}, Tuple{MeshArray{T, N, MT}, Vararg{Int64, N}}} where {T, MT, N}","page":"MeshArrays","title":"Base.getindex","text":"getindex(obj::MeshArray, inds...)\n\nReturn a subset of obj's data as specified by inds, where each inds may be, for example, an Int, an AbstractRange, or a Vector. \n\n\n\n\n\n","category":"method"},{"location":"lib/mesharrays/#Base.setindex!-Union{Tuple{N}, Tuple{MT}, Tuple{T}, Tuple{MeshArray{T, N, MT}, Any, Vararg{Int64, N}}} where {T, MT, N}","page":"MeshArrays","title":"Base.setindex!","text":"setindex!(obj::MeshArray, v, inds...)\nobj[inds...] = v\n\nStore values from array v within some subset of obj.data as specified by inds.\n\n\n\n\n\n","category":"method"},{"location":"lib/mesharrays/#Base.similar-Union{Tuple{S}, Tuple{N}, Tuple{MT}, Tuple{T}, Tuple{MeshArray{T, N, MT}, Type{S}}} where {T, MT, N, S}","page":"MeshArrays","title":"Base.similar","text":"Base.similar(obj::MeshArray{T,N,MT}, ::Type{S}) where {T,MT,N,S}\nBase.similar(obj::MeshArray{T,N,MT}) where {T,MT,N} = Base.similar(obj, T)\n\nReturn type:\n\nBase.similar(obj::MeshArray): Return a new MeshArray with the same meshes, and the uninitialized data of the same type as obj.data.\nBase.similar(obj::MeshArray, ::Type{S}): Return a new MeshArray with the same meshes, but with the uninitialized data of type S.\n\n\n\n\n\n","category":"method"},{"location":"lib/mesharrays/#Base.size-Tuple{GreenFunc.MeshArrays.AbstractMeshArray}","page":"MeshArrays","title":"Base.size","text":"size(obj::AbstractMeshArray)\n\nReturn a tuple containing the dimensions of obj.data (obj.dims).\n\n\n\n\n\n","category":"method"},{"location":"lib/mesharrays/#GreenFunc.MeshArrays._check-Tuple{MeshArray, MeshArray}","page":"MeshArrays","title":"GreenFunc.MeshArrays._check","text":"function _check(objL::MeshArray, objR::MeshArray)\n\nCheck if the Green's functions objL and objR are on the same meshes. Throw an AssertionError if any check is false.\n\n\n\n\n\n","category":"method"},{"location":"lib/mesharrays/#GreenFunc.MeshArrays.rank-Union{Tuple{Type{MeshArray{T, N, MT}}}, Tuple{N}, Tuple{MT}, Tuple{T}} where {T, MT, N}","page":"MeshArrays","title":"GreenFunc.MeshArrays.rank","text":"function rank(obj::MeshArray)\n\nReturn the dimension of obj.data (N).\n\n\n\n\n\n","category":"method"},{"location":"lib/greenfunc/#Green's-Functions","page":"GreenFunc","title":"Green's Functions","text":"","category":"section"},{"location":"lib/greenfunc/","page":"GreenFunc","title":"GreenFunc","text":"Modules = [GreenFunc]","category":"page"},{"location":"lib/greenfunc/#Base.:<<-Union{Tuple{MT2}, Tuple{MT1}, Tuple{N}, Tuple{T}, Tuple{MeshArray{T, N, MT1}, MeshArray{T, N, MT2}}} where {T, N, MT1, MT2}","page":"GreenFunc","title":"Base.:<<","text":"Base.:<<(objL::MeshArray, objR::MeshArray)\n\nDLR Fourier transform of functions on the first temporal grid (ImTime, ImFreq or DLRFreq). \n\nIf objL and objR have identical temporal grid, objL<<objR assign objR to objL.\nIf objL and objR have different temporal grid, one of them has to be in DLR space.\nIf objL is in DLR space, objL<<objR calculates the DLR spectral density of data in objR\nif objR is in DLR space, objL<<objR calculates the Green's function from the DLR spectral density in objR.\n\n\n\n\n\n","category":"method"},{"location":"lib/greenfunc/#GreenFunc.dlr_to_imtime-Union{Tuple{MeshArray{T, N, MT}}, Tuple{MT}, Tuple{N}, Tuple{T}, Tuple{MeshArray{T, N, MT}, Union{Nothing, AbstractVector}}} where {T, N, MT}","page":"GreenFunc","title":"GreenFunc.dlr_to_imtime","text":"function dlr_to_imfreq(mesharray[, tgrid; dim])\n\nTransform a Green's function in DLR to the imaginary-time domain.  #Arguements\n\n'mesharray': MeshArray in DLR space\ntgrid: The imaginary-time grid which the function transforms into. Default value is the imaginary-time grid from the DLRGrid from mesharray.mesh[dim].\ndim: The dimension of the temporal mesh. Default value is the first ImTime mesh.\n\n\n\n\n\n","category":"method"},{"location":"lib/greenfunc/#GreenFunc.imfreq_to_dlr-Union{Tuple{MeshArray{T, N, MT}}, Tuple{MT}, Tuple{N}, Tuple{T}} where {T, N, MT}","page":"GreenFunc","title":"GreenFunc.imfreq_to_dlr","text":"function imfreq_to_dlr(mesharray[; dim])\n\nCalculate the DLR sepctral density of a Matsubara-frequency Green's function. #Arguements\n\n'mesharray': MeshArray in the Matsubara-frequency domain.\ndim: The dimension of the mesh to be transformed. Default value is the first dimension with mesh type ImFreq.\n\n\n\n\n\n","category":"method"},{"location":"lib/greenfunc/#GreenFunc.imtime_to_dlr-Union{Tuple{MeshArray{T, N, MT}}, Tuple{MT}, Tuple{N}, Tuple{T}} where {T, N, MT}","page":"GreenFunc","title":"GreenFunc.imtime_to_dlr","text":"function imtime_to_dlr(mesharray[; dim])\n\nCalculate the DLR sepctral density of an imaginary-time Green's function.\n\n#Arguements\n\n'mesharray': MeshArray in the imaginary-time domain.\ndim: The dimension of the mesh to be transformed. Default value is the first dimension with mesh type ImTime.\n\n\n\n\n\n","category":"method"},{"location":"lib/greenfunc/#GreenFunc.to_dlr-Union{Tuple{MeshArray{T, N, MT}}, Tuple{MT}, Tuple{N}, Tuple{T}} where {T, N, MT}","page":"GreenFunc","title":"GreenFunc.to_dlr","text":"function to_dlr(mesharray[; dim])\n\nCalculate the DLR sepctral density of an imaginary-time or Matsubara-frequency Green's function.\n\n#Arguements\n\n'mesharray': MeshArray in the imaginary-time or the Matsubara-frequency domain.\ndim: The dimension of the mesh to be transformed. Default value is the first dimension with mesh type ImTime or ImFreq.\n\n\n\n\n\n","category":"method"},{"location":"lib/greenfunc/#GreenFunc.to_imfreq-Union{Tuple{MeshArray{T, N, MT}}, Tuple{MT}, Tuple{N}, Tuple{T}} where {T, N, MT}","page":"GreenFunc","title":"GreenFunc.to_imfreq","text":"function to_imfreq(mesharray[; dim])\n\nTransform a Green's function to the Matsubara-frequency domain.\n\n#Arguements\n\n'mesharray': MeshArray in the imaginary-time, the Matsubara-frequency or the DLR frequency domain.\ndim: The dimension of the mesh to be transformed. Default value is the first dimension with mesh type DLRFreq, ImTime or ImFreq.\n\n\n\n\n\n","category":"method"},{"location":"lib/greenfunc/#GreenFunc.to_imtime-Union{Tuple{MeshArray{T, N, MT}}, Tuple{MT}, Tuple{N}, Tuple{T}} where {T, N, MT}","page":"GreenFunc","title":"GreenFunc.to_imtime","text":"function to_imtime(mesharray[; dim])\n\nTransform a Green's function to the imaginary-time domain.\n\n#Arguements\n\n'mesharray': MeshArray in the imaginary-time, the Matsubara-frequency or the DLR frequency domain.\ndim: The dimension of the mesh to be transformed. Default value is the first dimension with mesh type DLRFreq, ImTime or ImFreq.\n\n\n\n\n\n","category":"method"},{"location":"lib/meshgrids/#MeshGrids","page":"MeshGrids","title":"MeshGrids","text":"","category":"section"},{"location":"lib/meshgrids/","page":"MeshGrids","title":"MeshGrids","text":"Modules = [GreenFunc.MeshGrids]","category":"page"},{"location":"lib/meshgrids/#GreenFunc.MeshGrids.DLRFreq","page":"MeshGrids","title":"GreenFunc.MeshGrids.DLRFreq","text":"function DLRFreq(β, isFermi::Bool=false;\n    dtype=Float64,\n    rtol=1e-12,\n    Euv=1000 / β,\n    sym=:none,\n    rebuild=false,\n    dlr::Union{DLRGrid,Nothing}=nothing\n)\n\nCreate a DLRFreq struct from parameters.\n\nArguments\n\nβ: inverse temperature.\nisFermi: the statistics for particles is fermionic or not. False by default.\ndtype: type of β and Euv.\nrtol: tolerance absolute error. By default, rtol = 1e-12.\nEuv: the UV energy scale of the spectral density. By default, Euv = 1000 / β.\nsymmetry: :ph for particle-hole symmetric, :pha for particle-hole symmetry, and :none for no symmetry. By default, sym = :none.\nrebuild: if no dlr is input, set false to load DLRGrid from the file; set true to recalculate the DLRGrid on the fly. By default, rebuild = false.\n\n\n\n\n\n","category":"type"},{"location":"lib/meshgrids/#GreenFunc.MeshGrids.DLRFreq-2","page":"MeshGrids","title":"GreenFunc.MeshGrids.DLRFreq","text":"struct DLRFreq{T<:Real} <: TemporalGrid{Int}\n\nDiscrete-Lehmann-representation grid for Green's functions. \n\nParameters\n\nT: type of the grid point, β and Euv.\n\nMembers\n\ndlr: built-in DLR grid.\ngrid: 1D grid of time axis, with locate, volume, and AbstractArray interface implemented. It should be grid of Int for ImFreq, and DLRGrid for DLRFreq.\nβ: inverse temperature.\nEuv:  the UV energy scale of the spectral density.\nrtol: tolerance absolute error.\nsymmetry: :ph for particle-hole symmetric, :pha for particle-hole symmetry, and :none for no symmetry. By default, sym = :none.\nisFermi: the statistics for particles. \n\n\n\n\n\n","category":"type"},{"location":"lib/meshgrids/#GreenFunc.MeshGrids.DLRFreq-Tuple{Lehmann.DLRGrid}","page":"MeshGrids","title":"GreenFunc.MeshGrids.DLRFreq","text":"function DLRFreq(dlr::DLRGrid)\n\nCreate a DLRFreq struct from DLRGrid.\n\nArguments\n\ndlr: 1D DLR grid.\n\n\n\n\n\n","category":"method"},{"location":"lib/meshgrids/#GreenFunc.MeshGrids.ImFreq","page":"MeshGrids","title":"GreenFunc.MeshGrids.ImFreq","text":"function ImFreq(β, isFermi::Bool=false;\n    dtype=Float64,\n    Euv=1000 / β,\n    rtol=1e-12,\n    symmetry=:none,\n    grid::Union{AbstractGrid,AbstractVector,Nothing}=nothing\n)\n\nCreate a ImFreq struct from parameters.\n\nArguments\n\nβ: inverse temperature.\nisFermi: the statistics for particles is fermionic or not. False by default.\ndtype: type of β and Euv. By default, dtype = Float64.\nEuv: the UV energy scale of the spectral density. By default, Euv = 1000 / β.\nsymmetry: :ph for particle-hole symmetric, :pha for particle-hole symmetry, and :none for no symmetry. By default, sym = :none.\ngrid: 1D Matsubara-frequency integer-valued grid as a AbstractVector or CompositeGrids.AbstractGrid. By default, a optimized grid built in DLR is used.\n\n\n\n\n\n","category":"type"},{"location":"lib/meshgrids/#GreenFunc.MeshGrids.ImFreq-2","page":"MeshGrids","title":"GreenFunc.MeshGrids.ImFreq","text":"struct ImFreq{T, G, R} <: TemporalGrid{Int}\n\nImaginary-frequency grid for Green's functions. \n\nParameters\n\nT<:Real: type of the grid point, β and Euv.\nG<:AbstractGrid{T}: type of 1D grid with T as the grid point type.\nR: type of the representation.\n\nMembers\n\ngrid: 1D grid of time axis, with locate, volume, and AbstractArray interface implemented. It should be grid of Int for ImFreq, and DLRGrid for DLRFreq.\nβ: inverse temperature.\nEuv:  the UV energy scale of the spectral density.\nisFermi: the statistics for particles is fermionic or not.\nsymmetry: :ph for particle-hole symmetric, :pha for particle-hole symmetry, and :none for no symmetry. By default, sym = :none.\nrtol: relative tolerance\nrepresentation: the representation of the Green's function.\n\n\n\n\n\n","category":"type"},{"location":"lib/meshgrids/#GreenFunc.MeshGrids.ImFreq-Tuple{Lehmann.DLRGrid}","page":"MeshGrids","title":"GreenFunc.MeshGrids.ImFreq","text":"function ImFreq(dlr::DLRGrid; dtype=Float64, grid::Union{AbstractGrid,AbstractVector}=SimpleG.Arbitrary{Int}(dlr.n))\n\nConstruct ImFreq from a DLRGrid, with a given grid. By default, grid is the Matsubara-frequency points from DLRGrid.\n\n\n\n\n\n","category":"method"},{"location":"lib/meshgrids/#GreenFunc.MeshGrids.ImTime","page":"MeshGrids","title":"GreenFunc.MeshGrids.ImTime","text":"function ImTime(β, isFermi::Bool=false;\n    dtype=Float64,\n    rtol=1e-12,\n    Euv=1000 / β,\n    symmetry=:none,\n    grid::Union{AbstractGrid,AbstractVector,Nothing}=nothing\n)\n\nCreate a ImTime struct.\n\nArguments\n\nβ: inverse temperature.\nisFermi: the statistics for particles is fermionic or not. False by default.\ndtype: type of the grid point. By default, dtype = Float64.\nEuv: the UV energy scale of the spectral density. By default, Euv = 1000 / β.\nsymmetry: :ph for particle-hole symmetric, :pha for particle-hole symmetry, and :none for no symmetry. By default, sym = :none.\ngrid: 1D time grid as a AbstractVector or CompositeGrids.AbstractGrid. By default, a optimized grid built in DLR is used.\n\n\n\n\n\n","category":"type"},{"location":"lib/meshgrids/#GreenFunc.MeshGrids.ImTime-2","page":"MeshGrids","title":"GreenFunc.MeshGrids.ImTime","text":"struct ImTime{T, G, R} <: TemporalGrid{T}\n\nTime grid for Green's functions.\n\nParameters\n\nT<:Real: type of the grid point, β and Euv.\nG<:AbstractGrid{T}: type of 1D grid with T as the grid point type.\nR: type of the representation.\n\nMembers\n\ngrid: 1D grid of time axis, with locate, volume, and AbstractArray interface implemented. It should be grid of Int for ImFreq, and DLRGrid for DLRFreq.\nβ: inverse temperature.\nEuv:  the UV energy scale of the spectral density.\nisFermi: the statistics for particles is fermionic or not.\nsymmetry: :ph for particle-hole symmetric, :pha for particle-hole symmetry, and :none for no symmetry. By default, sym = :none.\nrtol: relative tolerance\nrepresentation: the representation of the Green's function.\n\n\n\n\n\n","category":"type"},{"location":"lib/meshgrids/#GreenFunc.MeshGrids.ImTime-Tuple{Lehmann.DLRGrid}","page":"MeshGrids","title":"GreenFunc.MeshGrids.ImTime","text":"function ImTime(dlr::DLRGrid; dtype=Float64, grid::Union{AbstractGrid,AbstractVector}=SimpleG.Arbitrary{dtype}(dlr.τ))\n\nConstruct ImTime from a DLRGrid, with a given grid. By default, grid is the imaginary-time grid points from DLRGrid.\n\n\n\n\n\n","category":"method"},{"location":"lib/meshgrids/#GreenFunc.MeshGrids.MeshProduct","page":"MeshGrids","title":"GreenFunc.MeshGrids.MeshProduct","text":"The cartesian Mesh product:\n\n#Parameters:\n\n'MT': Type of meshes \n'N' : Number of meshes\n\n#Members:\n\n'meshes' : The list of Meshes in the MeshProduct\n'dims' : A tuple of the length of the mesh factors\n\n\n\n\n\n","category":"type"},{"location":"lib/meshgrids/#Base.getindex-Tuple{ImFreq, Int64}","page":"MeshGrids","title":"Base.getindex","text":"getindex(g::ImFreq, I::Int)\n\nEquivalent to g[I], get the real-valued Matsubara frequency of the Ith point in the grid.  For fermion, return (2g[I]+1)π/β, for boson, return 2g[I]*π/β.\n\nIf you need the integer-valued frequency, use g.grid[I] instead.\n\n\n\n\n\n","category":"method"},{"location":"lib/meshgrids/#Base.length-Tuple{MeshProduct}","page":"MeshGrids","title":"Base.length","text":"function Base.length(obj::MeshProduct)\n\nReturn the number of grids of the MeshProduct.\n\n\n\n\n\n","category":"method"},{"location":"lib/meshgrids/#Base.show-Tuple{IO, DLRFreq}","page":"MeshGrids","title":"Base.show","text":"show(io::IO, tg::DLRFreq)\n\nWrite a text representation of the DLR grid tg to the output stream io.\n\n\n\n\n\n","category":"method"},{"location":"lib/meshgrids/#Base.show-Tuple{IO, ImFreq}","page":"MeshGrids","title":"Base.show","text":"show(io::IO, tg::ImFreq)\n\nWrite a text representation of the Imaginary-frequency grid tg to the output stream io.\n\n\n\n\n\n","category":"method"},{"location":"lib/meshgrids/#Base.show-Tuple{IO, ImTime}","page":"MeshGrids","title":"Base.show","text":"show(io::IO, tg::ImTime)\n\nWrite a text representation of the Imaginary-time grid tg to the output stream io.\n\n\n\n\n\n","category":"method"},{"location":"lib/meshgrids/#Base.show-Tuple{IO, MeshProduct}","page":"MeshGrids","title":"Base.show","text":"function Base.show(io::IO, obj::MeshProduct)\n\nPrint the MeshProduct.\n\n\n\n\n\n","category":"method"},{"location":"lib/meshgrids/#Base.size-Tuple{MeshProduct, Int64}","page":"MeshGrids","title":"Base.size","text":"function Base.size(obj::MeshProduct, I::Int)\n\nReturn the length of the specifict Ith mesh factor of the MeshProduct.\n\n\n\n\n\n","category":"method"},{"location":"lib/meshgrids/#Base.size-Tuple{MeshProduct}","page":"MeshGrids","title":"Base.size","text":"function Base.size(obj::MeshProduct, I::Int)\n\nReturn the length of the specifict Ith mesh factor of the MeshProduct.\n\n\n\n\n\n","category":"method"},{"location":"lib/meshgrids/#GreenFunc.MeshGrids.index_to_linear-Union{Tuple{N}, Tuple{MT}, Tuple{MeshProduct{MT, N}, Vararg{Any}}} where {MT, N}","page":"MeshGrids","title":"GreenFunc.MeshGrids.index_to_linear","text":"function index_to_linear(obj::MeshProduct, index...)\n\nConvert a tuple of the indexes of each mesh to a single linear index of the MeshProduct.\n\nArgument:\n\n'obj': The MeshProduct object\n'index...': N indexes of the mesh factor, where N is the number of mesh factor\n\n\n\n\n\n","category":"method"},{"location":"lib/meshgrids/#GreenFunc.MeshGrids.linear_to_index-Union{Tuple{N}, Tuple{MT}, Tuple{MeshProduct{MT, N}, Int64}} where {MT, N}","page":"MeshGrids","title":"GreenFunc.MeshGrids.linear_to_index","text":"function linear_to_index(obj::MeshProduct, I::Int)\n\nConvert the single linear index of the MeshProduct to a tuple of indexes of each mesh. \n\nArgument:\n\n'obj': The MeshProduct object\n'I': The linear index of the MeshProduct \n\n\n\n\n\n","category":"method"},{"location":"lib/meshgrids/#GreenFunc.MeshGrids.locate-Tuple{ImFreq, Int64}","page":"MeshGrids","title":"GreenFunc.MeshGrids.locate","text":"locate(tg::ImFreq, n::Int)\nlocate(tg::ImFreq, ωn)\n\nFind the location in tg.grid for the Matsubara frequency ωn or the integer n.\n\n\n\n\n\n","category":"method"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = GreenFunc","category":"page"},{"location":"#GreenFunc","page":"Home","title":"GreenFunc","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for GreenFunc.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"lib/triqs/#Interface-to-Triqs","page":"Triqs","title":"Interface to Triqs","text":"","category":"section"},{"location":"lib/triqs/","page":"Triqs","title":"Triqs","text":"Modules = [GreenFunc.Triqs]","category":"page"},{"location":"lib/triqs/#GreenFunc.MeshArrays.MeshArray-Tuple{PythonCall.Py}","page":"Triqs","title":"GreenFunc.MeshArrays.MeshArray","text":"function MeshArray(objSrc::Py)\n\nConvert a Green's function object from triqs to a MeshArray.\n\n\n\n\n\n","category":"method"},{"location":"lib/triqs/#GreenFunc.Triqs.from_triqs-Tuple{PythonCall.Py}","page":"Triqs","title":"GreenFunc.Triqs.from_triqs","text":"function from_triqs(pyobj::Py)\n\nConvert a triqs object to a julia object. Currently support the following types:\n\nTriqs Mesh (ImTime, ImFreq, MeshProduct) -> MeshGrids (ImTime, ImFreq, MeshProduct)\nTriqs Green's function (Gf, GfImTime, GfImFreq)  -> MeshArray\nTriqs BlockGf -> Dict{String, MeshArray}\n\n\n\n\n\n","category":"method"}]
}
