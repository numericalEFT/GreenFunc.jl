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
mutable struct MeshArray{T,N,MT} <: AbstractMeshArray{T,N}
    #########   Mesh   ##############
    mesh::MT
    data::Array{T,N}
    dims::NTuple{N,Int}
end

"""
    function MeshArray(mesh...;
        innerstate::Union{AbstractVector{Int},Tuple{Vararg{Int}}}=(),
        data::Union{Nothing,AbstractArray}=nothing) where {T}
    
Create a Green struct. Its memeber `dims` is setted as the tuple consisting of the length of all meshes.

# Arguments
- `mesh`: meshes of Green's function. Mesh could be any iterable object, examples are vector, tuple, array, number, UnitRange (say, 1:5).
- `dtype`: data type of Green's function's value.
- `data`: the data of the Green's function. By default, `data` is constructed to an unintialized Array with the `dims` size containing elements of `dtype`.
"""
function MeshArray(mesh...;
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
    return MeshArray{dtype,N,typeof(mesh)}(mesh, data, dims)
end

########## Array Interface: https://docs.julialang.org/en/v1/manual/interfaces/#man-interface-array #############

"""
    size(obj::MeshArray)

Return a tuple containing the dimensions of `obj.data` (`obj.dims`).
"""
Base.size(obj::MeshArray) = obj.dims

"""
    eltype(obj::MeshArray)

Return the type of the elements contained in `obj.data`.
"""
Base.eltype(::Type{MeshArray{T,N,MT}}) where {T,N,MT} = T


"""
    getindex(obj::MeshArray, inds...)

Return a subset of `obj`'s data as specified by `inds`, where each `inds` may be, for example, an Int, an AbstractRange, or a Vector. 
"""
Base.getindex(obj::MeshArray{T,N,MT}, inds::Vararg{Int,N}) where {T,MT,N} = Base.getindex(obj.data, inds...)
# Base.getindex(obj::MeshArray, I::Int) = Base.getindex(obj.data, I)

"""
    setindex!(obj::MeshArray, v, inds...)
    obj[inds...] = v

Store values from array `v` within some subset of `obj.data` as specified by `inds`.
"""
Base.setindex!(obj::MeshArray{T,N,MT}, v, inds::Vararg{Int,N}) where {T,MT,N} = Base.setindex!(obj.data, v, inds...)
# Base.setindex!(obj::MeshArray, v, I::Int) = Base.setindex!(obj.data, v, I)




# IndexStyle(::Type{<:MeshArray}) = IndexCartesian() # by default, it is IndexCartesian

"""
    Base.similar(obj::MeshArray{T,N,MT}, ::Type{S}) where {T,MT,N,S}
    Base.similar(obj::MeshArray{T,N,MT}) where {T,MT,N} = Base.similar(obj, T)

# Return type:
- `Base.similar(obj::MeshArray)`: Return a new MeshArray with the same meshes, and the uninitialized data of the same type as `obj.data`.
- `Base.similar(obj::MeshArray, ::Type{S})`: Return a new MeshArray with the same meshes, but with the uninitialized data of type `S`.
"""
function Base.similar(obj::MeshArray{T,N,MT}, ::Type{S}) where {T,MT,N,S}
    return MeshArray(obj.mesh...; dtype=S, data=similar(obj.data, S))
end
Base.similar(obj::MeshArray{T,N,MT}) where {T,MT,N} = Base.similar(obj, T)
#By default, the following functions will all call Base.similar(obj::MeshArray, ::Type{S}, inds) as explained in https://docs.julialang.org/en/v1/manual/interfaces/#man-interface-array
#`Base.similar(obj::MeshArray, ::Type{S}, inds)`: Return a slice of obj.data.
#However, we don't want that since slice of GreeNew itself is not well defined with meshes.

################################ broadcast interface ###############################################
Base.BroadcastStyle(::Type{<:MeshArray}) = Broadcast.ArrayStyle{MeshArray}()

function Base.similar(bc::Base.Broadcast.Broadcasted{Broadcast.ArrayStyle{MeshArray}}, ::Type{ElType}) where {ElType}
    # println("get called")
    # Scan the inputs for the MeshArray:
    A = find_gf(bc)
    # Use other fields of A to create the output
    MeshArray(A.mesh, similar(Array{ElType}, axes(bc)), A.dims)
end

find_gf(bc::Broadcast.Broadcasted) = find_gf(bc.args)
find_gf(args::Tuple) = find_gf(find_gf(args[1]), Base.tail(args))
find_gf(x) = x
find_gf(::Tuple{}) = nothing
find_gf(a::MeshArray, rest) = a
find_gf(::Any, rest) = find_gf(rest)

function Base.copyto!(dest, bc::Base.Broadcast.Broadcasted{MeshArray{T,N,MT}}) where {T,MT,N}
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
# function Base.copyto!(dest::MeshArray{T, N, MT}, bc::Base.Broadcast.Broadcasted{Nothing}) where {T,MT,N}
#     _bcf = Base.Broadcast.flatten(bc)
#     bcf = Base.Broadcast.preprocess(dest, _bcf)
#     for I in CartesianIndices(dest)
#         dest[I] = bcf[I]
#     end
#     return dest
# end


# somehow, the following leads to stackoverflow due to some kind of infinite loop
# function Base.getproperty(obj::MeshArray{T,MT,N,Ninner}, sym::Symbol) where {T,MT,N,Ninner}
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
    show(io::IO, obj::MeshArray)

Write a text representation of the Green's function `obj` to the output stream `io`.
"""
function Base.show(io::IO, obj::MeshArray)
    print(io, "Green's function with dims = $(obj.dims) and total length = $(length(obj.data))\n"
              *
              "- Mesh: $(typeof(obj.mesh)) \n"
    )
end

"""
    function rank(obj::MeshArray{T,N,MT})

Return the dimension of `obj.data` (`N`).
"""
rank(::Type{MeshArray{T,N,MT}}) where {T,MT,N} = N


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
    function _check(objL::MeshArray, objR::MeshArray)

Check if the Green's functions `objL` and `objR` are on the same meshes. Throw an AssertionError if any check is false.
"""
function _check(objL::MeshArray, objR::MeshArray)
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
