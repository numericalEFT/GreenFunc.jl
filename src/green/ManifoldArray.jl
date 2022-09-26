# abstract type AbstractGreen{T,N,NINNER} <: AbstractArray{T,N} end

isiterable(::Type{T}) where {T} = hasmethod(iterate, (T,))
# isiterable(T) = hasmethod(iterate, (typeof(T),))

# """
#     mutable struct GreenDLR{T,Domain<:TimeDomain,TGT,MT,Ndata}

# General Green's function on a multi-dimensional mesh plus one in-built Discrete Lehmann Representation.

# # Parameters:
# - `T`: type of data
# - `Domain`: type of time domain, `Domain`<:`TimeDomain`.
# - `TGT`: type of time grid
# - `MT`: type of mesh
# - `N`: number of internal degrees of freedom
# - `Ndata`: rank of Green's function data, which always equals to N+2, 2 stands for the mesh and the extra dimension that has built-in DLR grid.

# # Members:
# - `DLR`: built-in DLR grid. Only one-dimensional DLR is available currently.
# - `tgrid` (TGT): the imaginary-time or Matsubara-frequency grid of dimension with built in DLR . If not provided by user, the optimized grid from DLR is used.
# - `mesh` (MT): the mesh is a direct product of grids of all other continuous degrees of freedom of Green's function, other than the one with DLR. The mesh has to support all standard Base functions of AbstractArray, plus the following two:
#     - locate(`mesh`, value): find the index of the closest grid point for given value;
#     - volume(`mesh`, index): find the volume of grid space near the point at griven index.
#     - volume(`mesh`, gridvalue): locate the corresponding index of a given value and than find the volume of grid space. 
# - `innerstate` (Tuple): innerstate saves the discrete inner dgrees of freedom of Green's function. 
# - `data` (Array{T,Ndata}): the data of the Green's function.
# """
# mutable struct ManifoldArray{T,MT,N} <: AbstractGreen{T,N}
mutable struct ManifoldArray{T,N,MT} <: AbstractArray{T,N}
    #########   Mesh   ##############
    mesh::MT
    data::Array{T,N}
    dims::NTuple{N,Int}
end

"""
    function ManifoldArray{T}(mesh...;
        innerstate::Union{AbstractVector{Int},Tuple{Vararg{Int}}}=(),
        data::Union{Nothing,AbstractArray}=nothing) where {T}
    
Create a Green struct. Its memeber `dims` is setted as the tuple consisting of the length of all meshes.

# Arguments
- `mesh`: meshes of Green's function. Mesh could be any iterable object, examples are vector, tuple, array, number, UnitRange (say, 1:5).
- `dtype`: data type of Green's function's value.
- `data`: the data of the Green's function. By default, `data` is constructed to an unintialized Array with the `dims` size containing elements of `dtype`.
"""
function ManifoldArray(mesh...;
    dtype=Float64,
    data::Union{Nothing,AbstractArray}=nothing)

    @assert all(x -> isiterable(typeof(x)), mesh) "all meshes should be iterable"

    N = length(mesh)
    dims = tuple([length(v) for v in mesh]...)
    if isnothing(data) == false
        # @assert length(size(data)) == N
        @assert size(data) == dims
    else
        data = Array{dtype,N}(undef, dims...)
    end
    return ManifoldArray{dtype,N,typeof(mesh)}(mesh, data, dims)
end

########## Array Interface: https://docs.julialang.org/en/v1/manual/interfaces/#man-interface-array #############

"""
    size(obj::ManifoldArray)

Return a tuple containing the dimensions of `obj.data` (`obj.dims`).
"""
Base.size(obj::ManifoldArray) = obj.dims

"""
    eltype(obj::ManifoldArray)

Return the type of the elements contained in `obj.data`.
"""
Base.eltype(::Type{ManifoldArray{T,N,MT}}) where {T,N,MT} = T


"""
    getindex(obj::ManifoldArray, inds...)

Return a subset of `obj`'s data as specified by `inds`, where each `inds` may be, for example, an Int, an AbstractRange, or a Vector. 
"""
Base.getindex(obj::ManifoldArray{T,N,MT}, inds::Vararg{Int,N}) where {T,MT,N} = Base.getindex(obj.data, inds...)
# Base.getindex(obj::ManifoldArray, I::Int) = Base.getindex(obj.data, I)

"""
    setindex!(obj::ManifoldArray, v, inds...)
    obj[inds...] = v

Store values from array `v` within some subset of `obj.data` as specified by `inds`.
"""
Base.setindex!(obj::ManifoldArray{T,N,MT}, v, inds::Vararg{Int,N}) where {T,MT,N} = Base.setindex!(obj.data, v, inds...)
# Base.setindex!(obj::ManifoldArray, v, I::Int) = Base.setindex!(obj.data, v, I)




# IndexStyle(::Type{<:ManifoldArray}) = IndexCartesian() # by default, it is IndexCartesian

"""
    Base.similar(obj::ManifoldArray{T,N,MT}, ::Type{S}) where {T,MT,N,S}
    Base.similar(obj::ManifoldArray{T,N,MT}) where {T,MT,N} = Base.similar(obj, T)

# Return type:
- `Base.similar(obj::ManifoldArray)`: Return a new ManifoldArray with the same meshes, and the uninitialized data of the same type as `obj.data`.
- `Base.similar(obj::ManifoldArray, ::Type{S})`: Return a new ManifoldArray with the same meshes, but with the uninitialized data of type `S`.
"""
function Base.similar(obj::ManifoldArray{T,N,MT}, ::Type{S}) where {T,MT,N,S}
    return ManifoldArray(obj.mesh...; dtype=S, data=similar(obj.data, S))
end
Base.similar(obj::ManifoldArray{T,N,MT}) where {T,MT,N} = Base.similar(obj, T)
#By default, the following functions will all call Base.similar(obj::ManifoldArray, ::Type{S}, inds) as explained in https://docs.julialang.org/en/v1/manual/interfaces/#man-interface-array
#`Base.similar(obj::ManifoldArray, ::Type{S}, inds)`: Return a slice of obj.data.
#However, we don't want that since slice of GreeNew itself is not well defined with meshes.

################################ broadcast interface ###############################################
Base.BroadcastStyle(::Type{<:ManifoldArray}) = Broadcast.ArrayStyle{ManifoldArray}()

function Base.similar(bc::Base.Broadcast.Broadcasted{Broadcast.ArrayStyle{ManifoldArray}}, ::Type{ElType}) where {ElType}
    # println("get called")
    # Scan the inputs for the ManifoldArray:
    A = find_gf(bc)
    # Use other fields of A to create the output
    ManifoldArray(A.mesh, similar(Array{ElType}, axes(bc)), A.dims)
end

find_gf(bc::Broadcast.Broadcasted) = find_gf(bc.args)
find_gf(args::Tuple) = find_gf(find_gf(args[1]), Base.tail(args))
find_gf(x) = x
find_gf(::Tuple{}) = nothing
find_gf(a::ManifoldArray, rest) = a
find_gf(::Any, rest) = find_gf(rest)

function Base.copyto!(dest, bc::Base.Broadcast.Broadcasted{ManifoldArray{T,N,MT}}) where {T,MT,N}
    # without this function, inplace operation like g1 .+= g2 will make a lot of allocations
    # Please refer to the following posts for more details:
    # 1. manual on the interface: https://docs.julialang.org/en/v1/manual/interfaces/#extending-in-place-broadcast-2
    # 2. see the post: https://discourse.julialang.org/t/help-implementing-copyto-for-broadcasting/51204/3
    # 3. example from DataFrames.jl: https://github.com/JuliaData/DataFrames.jl/blob/main/src/other/broadcasting.jl#L193

    ######## approach 2: use materialize ########
    bcf = Base.Broadcast.materialize(bc)
    for I in CartesianIndices(dest)
        dest[I] = bcf[I]
    end
    return dest
end

########### alternative approach ######################
# function Base.copyto!(dest::ManifoldArray{T, N, MT}, bc::Base.Broadcast.Broadcasted{Nothing}) where {T,MT,N}
#     _bcf = Base.Broadcast.flatten(bc)
#     bcf = Base.Broadcast.preprocess(dest, _bcf)
#     for I in CartesianIndices(dest)
#         dest[I] = bcf[I]
#     end
#     return dest
# end


# somehow, the following leads to stackoverflow due to some kind of infinite loop
# function Base.getproperty(obj::ManifoldArray{T,MT,N,Ninner}, sym::Symbol) where {T,MT,N,Ninner}
#     if sym === :N
#         return N
#     elseif sym === :Ninner
#         return Ninner
#     elseif sym === :dims
#         return obj.dims
#     else # fallback to getfield
#         return getfield(obj, sym)
#     end
# end

"""
    show(io::IO, obj::ManifoldArray)

Write a text representation of the Green's function `obj` to the output stream `io`.
"""
function Base.show(io::IO, obj::ManifoldArray)
    print(io, "Green's function with dims = $(obj.dims) and total length = $(length(obj.data))\n"
              *
              "- Mesh: $(typeof(obj.mesh)) \n"
    )
end

"""
    function rank(obj::ManifoldArray{T,N,MT})

Return the dimension of `obj.data` (`N`).
"""
rank(::Type{ManifoldArray{T,N,MT}}) where {T,MT,N} = N


#TODO:Following triqs design, we want the following two things:
# G<<init_function will initiate G with the given function.
# G<<g2 will copy g2 into G.


# """
#    def __lshift__(self, A):
#         A can be two things:
#           * G << any_init will init the GFBloc with the initializer
#           * G << g2 where g2 is a GFBloc will copy g2 into self
#         if isinstance(A, Gf):
#             if self is not A: # otherwise it is useless AND does not work !!
#                 assert self.mesh == A.mesh, "Green function meshes are not compatible:\n  %s\nand\n  %s" % (self.mesh, A.mesh)
#                 self.copy_from(A)
#         elif isinstance(A, lazy_expressions.LazyExpr): # A is a lazy_expression made of GF, scalars, descriptors
#             A2 = descriptors.convert_scalar_to_const(A)
#             def e_t (x):
#                 if not isinstance(x, descriptors.Base): return x
#                 tmp = self.copy()
#                 x(tmp)
#                 return tmp
#             self.copy_from (lazy_expressions.eval_expr_with_context(e_t, A2) )
#         elif isinstance(A, lazy_expressions.LazyExprTerminal): #e.g. g<< SemiCircular (...)
#             self << lazy_expressions.LazyExpr(A)
#         elif descriptors.is_scalar(A): #in the case it is a scalar ....
#             self << lazy_expressions.LazyExpr(A)
#         else:
#             raise NotImplemented
#                     return self
# """

"""
    function _check(objL::ManifoldArray, objR::ManifoldArray)

Check if the Green's functions `objL` and `objR` are on the same meshes. Throw an AssertionError if any check is false.
"""
function _check(objL::ManifoldArray, objR::ManifoldArray)
    # KUN: check --> __check
    # first:  check typeof(objL.tgrid)==typeof(objR.tgrid) 
    # second: check length(objL.tgrid)
    # third:  hasmethod(objL.tgrid, isequal) --> assert
    # @assert objL.innerstate == objR.innerstate "Green's function innerstates are not inconsistent: $(objL.innerstate) and $(objR.innerstate)"
    @assert typeof(objL.mesh) == typeof(objR.mesh) "Green's function meshes' types are inconsistent: $(typeof(objL.mesh)) and $(typeof(objR.mesh))"
    @assert objL.dims == objR.dims "Green's function dims are inconsistent: $(objL.dims) and $(objR.dims)"

    return true
    # @assert objL.tgrid == objR.tgrid "Green's function time grids are not compatible:\n $(objL.tgrid)\nand\n $(objR.tgrid)"
    # @assert objL.mesh == objR.mesh "Green's function meshes are not compatible:\n $(objL.mesh)\nand\n $(objR.mesh)"
end

 """
    <<(Obj::GreenDLR, objSrc::Py)
    Converts the green function from triqs to ManifoldArray.
"""
    
# """
#     <<(Obj::GreenDLR, objSrc::Expr)
#     Obj << objSrc

# Initiate the Green's function `Obj` with the given function expression `objSrc`.
# TODO: Add other behaviors including:
# 1. Convert other type of GreenFunc to GreenDLR
# 2. Assign a GreenDLR to another GreenDLR
# 3. Convert triqs Green's function to GreenDLR
# 4. First convert triqs object to triqs green's function, than to GreenDLR
# """
# function Base.:<<(Obj::GreenDLR, objSrc::Expr)
#     # init version of <<
#     # more general version needed
#     for (id, d) in enumerate(Obj)
#         inds = ind2sub_gen(size(Obj), id)
#         p, ωn, n, τ = NaN, NaN, NaN, NaN
#         G = d
#         β = Obj.β
#         if Obj.domain == ImFreq
#             n = Obj.tgrid[inds[3]]
#             if Obj.isFermi
#                 ωn = π * (2 * n + 1) / β
#             else
#                 ωn = π * 2 * n * β
#             end
#         elseif Obj.domain == ImTime
#             τ = tgrid[inds[3]]
#         end
#         p = Obj.mesh[inds[2]]
 
#         m = Dict(
#             :G => G,
#             :p => p,
#             :ωn => ωn,
#             :n => n,
#             :τ => τ,
#             :β => β
#         )
#         Obj[id] = DictParser.evalwithdict(objSrc, m)
#     end

#     return nothing
# end

"""
    Base.:<<(objL::ManifoldArray, objR::ManifoldArray)

DLR Fourier transform of functions that has exactly one TemporalGrid(ImTime, ImFreq or DLRFreq) among the meshes. 
If objL and objR have identical TemporalGrid, objL<<objR assign objR to objL.
If objL and objR have different TemporalGrid, one of them has to be in DLR space.
If objL is in DLR space, objL<<objR calculates the DLR spectral density of data in objR
if objR is in DLR space, objL<<objR calculates the corresponding data from the DLR spectral density in objR.
"""

function Base.:<<(objL::ManifoldArray, objR::ManifoldArray)
    # init version of <<
    # more general version needed
    axes=[]
    for (mi,mesh) in enumerate(objL.mesh) 
        if(typeof(mesh) <: TemporalGrid)
            append!(axes,mi)
        else
            #TODO:add hashtable for mesh to support == operation 
            @assert typeof(mesh) == typeof(objR.mesh[mi]) "Meshes not involved in Fourier transform have to be identical" #should assert mesh == objR.mesh[mi] when == is defined
        end
    end
    @assert length(axes)==1 "Only one temporal mesh with built in DLR grid is allowed"
    typeL = typeof(objL.mesh[axes[1]])
    typeR = typeof(objR.mesh[axes[1]])
    meshL = objL.mesh[axes[1]]
    meshR = objR.mesh[axes[1]]

    if typeL==typeR
        error("Assignment still work in progress") #need == operation for mesh
        # if meshL == meshR
        #     objL.data=deepcopy(objR.data)
        # else
        #     error("Green's function can only be assigned to another with the same mesh")
        # end
    elseif typeL <: MeshGrids.DLRFreq
        if typeR <: MeshGrids.ImFreq
            objL.data = matfreq2dlr(meshL.dlr, objR.data, meshR.grid; axis=axes[1])
        elseif typeR <: MeshGrids.ImTime
            objL.data = tau2dlr(meshL.dlr, objR.data, meshR.grid; axis=axes[1])
        end
    elseif typeR <: MeshGrids.DLRFreq
        if typeL <: MeshGrids.ImFreq
            objL.data = dlr2matfreq(meshR.dlr, objR.data, meshL.grid; axis=axes[1])
        elseif typeL <: MeshGrids.ImTime
            objL.data = dlr2tau(meshR.dlr, objR.data, meshL.grid; axis=axes[1])
        end
    else
        error("One of the Grren's function has to be in DLRfreq space to do Fourier transform")
    end
end


"""
    function dlr_to_imtime(obj::ManifoldArray; kwargs...)

Convert sepctral density in DLR space to imaginary time space for function with exactly one DLRFreq grid among meshes.
#Arguements
#- 'obj': Function in DLR space
#- 'targetGrid': The imaginary time grid which the function transforms into. Default value is the imaginary time space from the DLR grid in obj.
"""
function dlr_to_imtime(obj::ManifoldArray; kwargs...)
    # init version of <<
    # more general version needed
    axes=[]
    for (mi,mesh) in enumerate(obj.mesh) 
        if(typeof(mesh) <: TemporalGrid)
            @assert typeof(mesh)<:MeshGrids.DLRFreq "Green's function has to be in DLRFreq space"
            append!(axes,mi)
        end
    end
    @assert length(axes)==1 "Only one temporal mesh with built in DLR grid is allowed"
    mesh = obj.mesh[axes[1]]
    if :targetGrid in keys(kwargs)
        tgrid=kwargs[:targetGrid]
        @assert typeof(tgrid)==MeshGrids.ImTime "Target grid has to be ImTime type"
    else
        tgrid=MeshGrids.ImTime(mesh.β,mesh.statistics)
    end

    mesh_copy=Tuple( i==axes[1] ? tgrid : obj.mesh[i] for i in 1:length(obj.dims)  )
    data = dlr2tau(mesh.dlr, obj.data, tgrid.grid; axis=axes[1])
    datatype = typeof(data[1])
    return ManifoldArray(mesh_copy...; dtype = datatype  ,data=data)
end

"""
    function dlr_to_imfreq(obj::ManifoldArray; kwargs...)

Convert sepctral density in DLR space to matsubara frequency space for function with exactly one DLRFreq grid among meshes.
#Arguements
#- 'obj': Function in DLR space
#- 'targetGrid': The matsubara frequency grid which the function transforms into. Default value is the imaginary time space from the DLR grid in obj.
"""

function dlr_to_imfreq(obj::ManifoldArray; kwargs...)
    # init version of <<
    # more general version needed
    axes=[]
    for (mi,mesh) in enumerate(obj.mesh)
        if(typeof(mesh) <: TemporalGrid)
            @assert typeof(mesh)<:MeshGrids.DLRFreq "Green's function has to be in DLRFreq space"
            append!(axes,mi)
        end
    end
    @assert length(axes)==1 "Only one temporal mesh with built in DLR grid is allowed"
    mesh = obj.mesh[axes[1]]
    if :targetGrid in keys(kwargs)
        tgrid=kwargs[:targetGrid]
        @assert typeof(tgrid)==MeshGrids.ImFreq "Target grid has to be MeshGrids.ImFreq type"
    else
        tgrid=MeshGrids.ImFreq(mesh.β,mesh.statistics)
    end
    mesh_copy=Tuple( i==axes[1] ? tgrid : obj.mesh[i] for i in 1:length(obj.dims)  )
    data = dlr2matfreq(mesh.dlr, obj.data, tgrid.grid; axis=axes[1])
    datatype = typeof(data[1])
    return ManifoldArray(mesh_copy...; dtype = datatype  ,data=data)
end

"""
    function to_dlr(obj::ManifoldArray; kwargs...)

Calculate the sepctral density of a function in imaginary time/matsubara frequency space, with exactly one DLRFreq grid among meshes.
#Arguements
#- 'obj': Function in imaginary time/matsubara frequency space
#- 'targetGrid': The DLRFreq grid which the function transforms into. Default value is a DLRFreq gridwith the same beta and statistics of the temporal grid in obj, with default rtol=1e-12 and Euv = 1000/beta. 
"""

function to_dlr(obj::ManifoldArray; kwargs...)
    # init version of <<
    # more general version needed
    axes=[]
    for (mi,mesh) in enumerate(obj.mesh)
        if(typeof(mesh) <: TemporalGrid)
            @assert typeof(mesh)<:MeshGrids.ImTime||typeof(mesh)<:MeshGrids.ImFreq "Green's function has to be in ImFreq or ImTime space"
            append!(axes,mi)
        end
    end
    @assert length(axes)==1 "Only one temporal mesh with built in DLR grid is allowed"
    mesh = obj.mesh[axes[1]]
    if :targetGrid in keys(kwargs)
        tgrid=kwargs[:targetGrid]
        @assert typeof(tgrid)==MeshGrids.DLRFreq "Target grid has to be DLRFreq type"
    else
        tgrid=MeshGrids.DLRFreq(mesh.β,mesh.statistics)
    end
    mesh_copy=Tuple( i==axes[1] ? tgrid : obj.mesh[i] for i in 1:length(obj.dims)  )
    if typeof(mesh)<:MeshGrids.ImFreq
        data = matfreq2dlr(tgrid.dlr, obj.data, mesh.grid; axis=axes[1])
    else typeof(mesh)<:MeshGrids.ImTime
        data = tau2dlr(tgrid.dlr, obj.data, mesh.grid; axis=axes[1])
    end
    datatype = typeof(data[1])
    return ManifoldArray(mesh_copy...; dtype = datatype  ,data=data)
end




# Return the single-particle density matrix of the Green's function `obj`.
# """
# function density(obj::ManifoldArray; kwargs...)
#     G_ins = toTau(obj, [obj.β,]).data .* (-1)
#     return selectdim(G_ins, ndims(G_ins), 1)
# end

#rMatrix transform of the target space of a matrix valued Greens function.
#Sets the current Greens function :math:`g_{ab}` to the matrix transform of :math:`G_{cd}`
#using the left and right transform matrices :math:`L_{ac}` and :math:`R_{db}`.
#.. math::
#g_{ab} = \sum_{cd} L_{ac} G_{cd} R_{db}


# """
#     def from_L_G_R(self, L, G, R):
#         Parameters
#         ----------
#         L : (a, c) ndarray
#             Left side transform matrix.
#         G : Gf matrix valued target_shape == (c, d)
#             Greens function to transform.
#         R : (d, b) ndarray
#             Right side transform matrix.
#         Notes
#         -----
#         Only implemented for Greens functions with a single mesh.


#         assert self.rank == 1, "Only implemented for Greens functions with one mesh"
#         assert self.target_rank == 2, "Matrix transform only valid for matrix valued Greens functions"

#         assert len(L.shape) == 2, "L needs to be two dimensional"
#         assert len(R.shape) == 2, "R needs to be two dimensional"

#         assert L.shape[1] == G.target_shape[0], "Dimension mismatch between L and G"
#         assert R.shape[0] == G.target_shape[1], "Dimension mismatch between G and R"

#         assert L.shape[0] == self.target_shape[0], "Dimension mismatch between L and self"
#         assert R.shape[1] == self.target_shape[1], "Dimension mismatch between R and self"

#         if not L.strides == sorted(L.strides):
#             L = L.copy(order='C')

#         if not R.strides == sorted(R.strides):
#             R = R.copy(order='C')

#         wrapped_aux.set_from_gf_data_mul_LR(self.data, L, G.data, R)
# """
# function from_L_G_R(self,L,G::ManifoldArray,R)
#     return 1
# end
